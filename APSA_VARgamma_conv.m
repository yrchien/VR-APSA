function    [outputVector,...
    errorVector,...
    coefficientVector,gamma] =   APSA_VARgamma_conv(desired,input,S,rho)

%   Affine_projection.m
%       Implements the Complex Affine-Projection algorithm for COMPLEX valued data.
%       (Algorithm 4.6 - book: Adaptive Filtering: Algorithms and Practical
%                                                       Implementation, Diniz)
%
%   Syntax:
%       [outputVector,errorVector,coefficientVector] =
%       Affine_projectionchina(desired,input,S,gamma)
%
%   Input Arguments:
%       . desired   : Desired signal.                               (ROW vector)
%       . input     : Signal fed into the adaptive filter.          (ROW vector)
%       . S         : Structure with the following fields
%           - step                  : Convergence (relaxation) factor.
%           - filterOrderNo         : Order of the FIR filter.
%           - initialCoefficients   : Initial filter coefficients.  (COLUMN vector)
%           -
%           - memoryLength          : Reuse data factor.
%                                     (referred as L in the textbook)

% gamma                 : Regularization factor.
%                                     (small positive constant to avoid singularity)


%   Output Arguments:
%       . outputVector      :   Store the estimated output of each iteration.   (COLUMN vector)
%       . errorVector       :   Store the error for each iteration.             (COLUMN vector)
%       . coefficientVector :   Store the estimated coefficients for each iteration.
%                               (Coefficients at one iteration are COLUMN vector)



%   Some Variables and Definitions:
%       . prefixedInput         :   Input is prefixed by nCoefficients -1 zeros.
%                                   (The prefix led to a more regular source code)
%
%       . regressor             :   Auxiliar variable. Store the piece of the
%                                   prefixedInput that will be multiplied by the
%                                   current set of coefficients.
%                                   (regressor is a COLUMN vector)
%
%       . nCoefficients         :   FIR filter number of coefficients.
%
%       . nIterations           :   Number of iterations.


%User defined parameter
gamma_min=1e-4;

%   Initialization Procedure
nCoefficients       =   S.filterOrderNo+1;
nIterations         =   length(desired);


%   Pre Allocations

errorVectorAp       =   zeros((S.memoryLength+1)    ,nIterations);
outputVectorAp      =   zeros((S.memoryLength+1)    ,nIterations);
coefficientVector       =   zeros(nCoefficients         ,(nIterations+1));
regressor               =   zeros(nCoefficients         ,(S.memoryLength+1));
gamma=zeros(1,(nIterations+1));
gamma(1)=1e-4;

%   Initial State Weight Vector
coefficientVector(:,1)  =   S.initialCoefficients;

%   Improve source code regularity


dvect=zeros((S.memoryLength+1),1);
xvect=zeros(nCoefficients,1);
xvectap =   zeros( (S.memoryLength+1),nCoefficients  );
old_regressor =   zeros( (S.memoryLength+1),nCoefficients  );

%   Body
for it = 1:nIterations,

    xvect=[input(it);xvect(1:(nCoefficients-1))];

    xvectap=[xvect';xvectap(1:S.memoryLength,:)];

    regressor=xvectap';

    dvect=[desired(it);dvect(1:S.memoryLength,:)];

    outputVectorAp(:,it)        =   (regressor')*coefficientVector(:,it);


    errorVectorAp(:,it)   =   (dvect)-outputVectorAp(:,it);

    %gamma(it+1)=gamma(it)-rho*sign(errorVectorAp(1,it).'*(regressor(:,1)')*(regressor*sign(errorVectorAp(:,it))));  %gamma
    if it>1
        Temp=old_regressor*(errorVectorAp(:,it-1));
        gamma(it+1)=gamma(it)-rho*sign(errorVectorAp(1,it).'*(regressor(:,1)')*(old_regressor*sign(errorVectorAp(:,it-1)))); %sign括號位置不同 效果差
    else
        gamma(it+1)=gamma(it);
    end

    if gamma(it+1)<gamma_min
        gamma(it+1)=gamma_min;
    end


    coefficientVector(:,it+1)       =   coefficientVector(:,it)+S.step*regressor*sign(errorVectorAp(:,it))/(sqrt(sign(errorVectorAp(:,it)')*(regressor')*regressor*(sign(errorVectorAp(:,it)))+gamma(it+1)));%有改gamma

    old_regressor= regressor;


end

outputVector            =   outputVectorAp(1,:)';
errorVector             =   errorVectorAp(1,:)';

%   EOF
