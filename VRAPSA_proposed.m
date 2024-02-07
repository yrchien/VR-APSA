function    [outputVector,errorVector,coefficientVector,gamma,delta_gamma] = VRAPSA_proposed(desired,input,S,rho)

% Initialization Procedure
nCoefficients=S.filterOrderNo+1;
nIterations=length(desired);

%   Pre-Allocations
errorVectorAp=zeros((S.memoryLength+1),nIterations); %P x iter
outputVectorAp=zeros((S.memoryLength+1),nIterations); % P x iter
coefficientVector=zeros(nCoefficients,(nIterations+1)); % L x iter
dvect=zeros((S.memoryLength+1),1); %P x 1
xvect=zeros(nCoefficients,1); %Lx1
xvectap=zeros((S.memoryLength+1),nCoefficients); %PxL
old_regressor=zeros((S.memoryLength+1),nCoefficients); %PxL

gamma=zeros(1,(nIterations+1));
delta_gamma=zeros(1,(nIterations+1));

%   Initial State Weight Vector
coefficientVector(:,1)=S.initialCoefficients;
gamma(1)=1e-4;
gamma_min=gamma(1);

%   Body
for it = 1:nIterations
    xvect=[input(it);xvect(1:end-1)];%Lx1
    xvectap=[xvect';xvectap(1:S.memoryLength,:)];%P x L
    regressor=xvectap'; % L x P
    dvect=[desired(it);dvect(1:end-1)];%P x 1
    outputVectorAp(:,it)=xvectap*coefficientVector(:,it); %P x 1
    errorVectorAp(:,it)=dvect-outputVectorAp(:,it);

    if it>1
        delta_gamma(it)=-rho*(regressor*sign(errorVectorAp(:,it))).'*old_regressor*sign(errorVectorAp(:,it-1));
        gamma(it+1)=gamma(it)+delta_gamma(it);  %gamma
    else
        gamma(it+1)=gamma(it);
        delta_gamma(it)=0;
    end


    %Set minimal value
    if gamma(it+1)<gamma_min
        gamma(it+1)=gamma_min;
    end

    % Coeff. Updating
    xs_vec=regressor*sign(errorVectorAp(:,it)); % L x 1
    coefficientVector(:,it+1)=coefficientVector(:,it)+S.step*xs_vec/(sqrt(xs_vec.'*xs_vec+gamma(it+1)));
    old_regressor= regressor;
end
outputVector=outputVectorAp(1,:)';
errorVector=errorVectorAp(1,:)';