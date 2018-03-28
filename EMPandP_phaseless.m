function [ r_hat ] = EMPandP_phaseless( A_gpu, At_gpu, y_gpu ,sigma_w,n,beta,AhA_diag,rho,iters,v_GS_init)
%EMPANDP Perform expectation maximization plug and play admm to recover the
% target object's reflectivity using an implicit (denoiser-imposed) prior
% from high SNR phaseless measurements
%   Forward model: y=|Ag+w_1|+w_2\approx exp(i\theta)\circ A g + w_2. p(g|r)=CN(0,diag(r))
%   Inputs: 
%   Outputs:
% See "Optically coherent image formation and denoising using a plug and play inversion framework"

[m,L]=size(y_gpu);
if L~=1
    error('Code currently only supports a single observation');
end

%%Initialization
if isempty(v_GS_init)
    iAhA=1./AhA_diag;
    Apinv_func=@(x) iAhA(:).*At_gpu(x);%Assumes AtA_inv=I
    v_GS_init=gather(GS(y_gpu, A_gpu, Apinv_func,100));
end
v=abs(v_GS_init).^2;
r=v;
u=zeros(n,1);
sigma_lambda=1/2*std(v);
theta=gather(angle(A_gpu(gpuArray(v_GS_init))));

%% ADMM Main Loop
for i=1:iters
    display(num2str(i));
    r_tilde=v-1/rho*u;
    [r,theta]=fidelity_proxmap(r_tilde,r,A_gpu, At_gpu,y_gpu,sigma_w,sigma_lambda,m,n,L,theta,AhA_diag);
    v_tilde=r+1/rho*u;
%     sigma_n=sqrt(beta*sigma_lambda^2);
    sigma_n=beta;%Currently I am just directly controlling the amount of smoothing
    v_tilde2=v_tilde.*(v_tilde>0);
%     v_tilde2(v_tilde2>1)=1;
    v=Denoise_BM3D(reshape(v_tilde2,[sqrt(n),sqrt(n)]),sigma_n);
    v=v(:);
    figure(1);
    subplot(2,2,1);imshow(abs(reshape(r_tilde,[sqrt(n),sqrt(n)])),[]);title('$\tilde{r}$');
    subplot(2,2,2);imshow(abs(reshape(r,[sqrt(n),sqrt(n)])),[]);title('r');
    subplot(2,2,3);imshow(abs(reshape(v_tilde,[sqrt(n),sqrt(n)])),[]);title('$\tilde{v}$');
    subplot(2,2,4);imshow(abs(reshape(v,[sqrt(n),sqrt(n)])),[]);title('v');
    u=u+rho*(r-v);
end
r_hat=v;

end
function [r_k,theta]=fidelity_proxmap(r_tilde,r_k,A, At,y_gpu,sigma_w,sigma_lambda,m,n,L,theta,AhA_diag)
%Perform MAP estimation over both the phase and the intensity. Use AM
    for iter=1:5
        
%         [Aty,~]=At_y_withPhase(At,y,theta);
        Aty=gather(At(exp(-1i*theta).*y_gpu)); 
%         C=diag(sparse(1./(1/sigma_w^2*L*AhA_diag.*ones(n,1)+1./(r_k+eps))));%Equivalent to above when A'A=AhA_diag*L*I. L=1 for my scenarios
        C_diag=1./(1/sigma_w^2*L*AhA_diag.*ones(n,1)+1./(r_k+eps));%Equivalent to above when A'A=AhA_diag*L*I. L=1 for my scenarios 
%         mu=C*1/sigma_w^2*Aty;
        mu=C_diag.*1/sigma_w^2.*Aty;
        
        %ML Estimate of phase shift
        z=A(mu);
        theta=gather(-angle(y_gpu.*z));
        
        
        %MAP Estimate of r
    %     C=inv(1/sigma_w^2*A'*A+diag(1./(r_k+eps)));
%         C=diag(sparse(1./(1/sigma_w^2.*ones(n,1)+1./(r_k+eps))));%Equivalent to above when A'A=I
%         C=diag(sparse(1./(1/sigma_w^2*L.*ones(n,1)+1./(r_k+eps))));%Equivalent to above when A'A=L*I
        Aty=gather(At(exp(-1i*theta).*y_gpu));
%         C=diag(sparse(1./(1/sigma_w^2*L*AhA_diag.*ones(n,1)+1./(r_k+eps))));%Equivalent to above when A'A=AhA_diag*L*I. L=1 for my scenarios 
        C_diag=1./(1/sigma_w^2*L*AhA_diag.*ones(n,1)+1./(r_k+eps));%Equivalent to above when A'A=AhA_diag*L*I. L=1 for my scenarios 
        
%         mu=C*1/sigma_w^2*Aty;
        mu=C_diag.*1/sigma_w^2.*Aty;
        

%         %Although roots supports gpuArray input arguments, it is not yet
%         %supported within arrayfun. https://www.mathworks.com/matlabcentral/answers/323056-functions-arrayfun-and-roots-with-gpu-computing
%         t1=tic;
%         alpha_array=gpuArray(zeros(4,n));
%         alpha_array(1,:)=1/sigma_lambda^2*ones(1,n);
%         alpha_array(2,:)=-r_tilde(:)'./sigma_lambda^2;
%         alpha_array(3,:)=ones(1,n);
%         alpha_array(4,:)=-(C_diag(:)'+abs(mu(:)').^2);
%         temp=arrayfun(@(x) roots(x),alpha_array);
%         temp_diffs=arrayfun(@(x) abs(x-real(x)),temp);
%         temp_i=arrayfun(@(x) x(x==min(x)),temp_diffs);
%         r_k=arrayfun(@(x) real(x(1)),temp_i);
%         toc(t1);

        t2=tic;
        parfor i=1:n
            alpha_1_i=1/sigma_lambda^2;
            alpha_2_i=-r_tilde(i)/sigma_lambda^2;
%             alpha_3_i=L;%This treats the signal as L indendent observations
            alpha_3_i=1;
%             alpha_4_i=-(L*C(i,i)+sum(abs(mu_ls(i,:)).^2,2));%This treats the signal as L indendent observations
%             alpha_4_i=-(C(i,i)+abs(mu(i)).^2);%This treats the signal as one observation, where the measurement matrix accounts for the phase shifts
            alpha_4_i=-(C_diag(i)+abs(mu(i)).^2);%This treats the signal as one observation, where the measurement matrix accounts for the phase shifts
            %alpha_1^3*r_i+alpha_2^2*r_i+alpha_3^1*r_i+alpha_4=0;
            alpha_i=[alpha_1_i,alpha_2_i,alpha_3_i,alpha_4_i];
            temp=roots(alpha_i);%Select among the roots the real valued one.
            temp_diffs=abs(temp-real(temp));
            temp_i=temp(temp_diffs==min(temp_diffs));
            r_k(i)=real(temp_i(1));
        end
        toc(t2)
        1;
    end
end
function out=Denoise_BM3D(x,sigma)
%     sigma=sigma*255;
    [~,out]=BM3D(1,x,sigma);
    out=255*out;
end
function [z,zs]=At_y_withPhase(At,y,phi)
    Aty=At(y);
    [n,L]=size(Aty);
    zs=zeros(n,L);
    for l=1:L
        zs(:,l)=exp(-1i*phi(l))*Aty(:,l);
    end
    z=sum(zs,2);
end
function [z,zs]=At_y_withPhase_zern(At,y,theta,psij)
    phase=exp(-1i*psij*theta);
    Aty=At(y);
    [n,L]=size(Aty);
    zs=zeros(n,L);
    for l=1:L
        zs(:,l)=phase(:,l).*Aty(:,l);
    end
    z=sum(zs,2);
end


