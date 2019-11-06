%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BELLA.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ x,f,out ] = BELLA(x0,func,kernel,subprob,options) 
% BELLA is a locally superlinearly convergent Bregman forward-backward  
% splitting method for solving nonsmooth nonconvex composite optimization
% 
%                    min f(x) + g(x) 
%      where x=(x_1,...,x_N) and
%            f is relatively smooth;
%            g is proper and lsc. 
%
% INPUT:
%
% func                 % function handle for the objective function
% kernel               % kernel of Bregman distance
% subprob              % function handle for associated subproblems
% x0                   % initial point
% options              % structure including the parameteres of scheme
%
%   .gamma0            % step-sizes for proximal operators 
%   .MaxNumIter        % maximum number of iterations
%   .MaxNumFunEval     % maximum number of function evaluations
%   .MaxNumGradEval    % maximum number of gradient evaluations
%   .TimeLimit         % maximum running time
%   .Stopping_Crit     % stopping criterion
%
%                      % 1 : stop if MaxNumIter is reached (default)
%                      % 2 : stop if MaxNumFunEval is reached
%                      % 3 : stop if MaxNumGradEval is reached
%                      % 4 : stop if TimeLimit is reached
%
% OUTPUT:
%
% x                    % the best approximation of the optimizer
% f                    % the best approximation of the optimum
% out                  % structure including more output information
%
%   .T                 % running time
%   .Niter             % total number of iterations
%   .Nfunc             % total number of function evaluations
%   .Ngrad             % total number of gradient evaluations
%   .F                 % array including all function values             
%   .Status            % reason of termination

% REFERENCE: 
%
% [1] Bregman forward-backward splitting for nonconvex composite optimization:
% superlinear convergence to nonisolated critical points, 
% Submitted,(2019)
%           
% WRITTEN BY: 
%
% Masoud Ahookhosh
% Department of Electrical Engineering(ESAT-STADIUS), KU Leuven, Belgium
%
% Yang Lu
% Department of Mathematics (Section of Statistics), KU Leuven, Belgium
%
% LAST UPDATE: 
%
% November 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ x,f,out ] = BELLA( x0,func,kernel,subprob,options )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Initializing and setting the parameters %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long;

% ================ Error messages for input and output =================
if nargin > 5
    error('The number of input arguments is more than what is needed');
elseif nargin < 5
    error('The number of input arguments is not enough');
end

if isempty(func)
    error('BPALM needs the function handle func to be defined');
elseif ~isa(func,'function_handle')
    error('func should be a function handle');
end

if isempty(subprob)
    error('BPALM needs the function handle subprob to be defined');
elseif ~isa(subprob,'function_handle')
    error('subprob should be a function handle');
end

% =================== initializing the parameters ======================
% ===== user has requested viewing the default values of "options" =====
[MaxNumIter,MaxNumFunEval,MaxNumGradEval,TimeLimit, ...
    flag_time,Stopping_Crit] = Initialization(options);

if flag_time == 1
    Time(1) = 0;
end

N        = length(x0);
gammak   = options.gamma0;
xk       = x0;
Niter    = 1;
fxk      = func(xk,0);
F        = fxk;
Nfunc    = 1;
Ngrad    = 0;
StopFlag = 0;

outx=[]; % Hien added

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Main body of BELLA.m %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T0 = tic;

% ======================= start of the main loop =======================
while ~StopFlag
    xk1 = xk;
        % solving subproblem %
        optsub.block   = i;
        gammaki        = gammak(i);
        optsub.gammaki = gammaki;
        [fxk1,gfki]    = func(xk1,i);
        Nfunc          = Nfunc+1;
        Ngrad          = Ngrad+1;
        [hxk1,ghki]    = kernel(xk1,i);
        optsub.gfki    = gfki;
        optsub.ghki    = ghki;
        xk1i           = subprob(xk1,optsub); 
        % updating the sequence %
        xk1(i)         = {xk1i};
        Niter          = Niter+1
        F              = [F,fxk1];
    xk = xk1;
    outx=[outx; xk1]; % Hien added
  
    % ================== checking stopping criteria ====================
    Time = toc(T0); 
    [StopFlag, Status] = StopCriterion(Niter,Nfunc,Ngrad, ...
    Time,MaxNumIter,MaxNumFunEval,MaxNumGradEval,TimeLimit, ...
    Stopping_Crit);   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Status
x          = xk;
f          = func(x,0);
T          = Time;
out.T      = T;
out.F      = F';
out.Niter  = Niter;
out.Nfunc  = Nfunc;
out.Ngrad  = Ngrad;
out.Status = Status;

out.x=outx; % Hien added

if flag_time == 1
    out.Time = Time;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% End of BELLA.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
