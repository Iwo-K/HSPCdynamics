% simulate_pd_fv.m is the matlab interface to the cvodes mex
%   which simulates the ordinary differential equation and respective
%   sensitivities according to user specifications.
%   this routine was generated using AMICI commit # in branch unknown branch in repo unknown repository.
%
% USAGE:
% ======
% [...] = simulate_pd_fv(tout,theta)
% [...] = simulate_pd_fv(tout,theta,kappa,data,options)
% [status,tout,x,y,sx,sy] = simulate_pd_fv(...)
%
% INPUTS:
% =======
% tout ... 1 dimensional vector of timepoints at which a solution to the ODE is desired
% theta ... 1 dimensional parameter vector of parameters for which sensitivities are desired.
%           this corresponds to the specification in model.sym.p
% kappa ... 1 dimensional parameter vector of parameters for which sensitivities are not desired.
%           this corresponds to the specification in model.sym.k
% data ... struct containing the following fields:
%     Y ... 2 dimensional matrix containing data.
%           columns must correspond to observables and rows to time-points
%     Sigma_Y ... 2 dimensional matrix containing standard deviation of data.
%           columns must correspond to observables and rows to time-points
%     T ... (optional) 2 dimensional matrix containing events.
%           columns must correspond to event-types and rows to possible event-times
%     Sigma_T ... (optional) 2 dimensional matrix containing standard deviation of events.
%           columns must correspond to event-types and rows to possible event-times
% options ... additional options to pass to the cvodes solver. Refer to the cvodes guide for more documentation.
%    .atol ... absolute tolerance for the solver. default is specified in the user-provided syms function.
%        default value is 1e-16
%    .rtol ... relative tolerance for the solver. default is specified in the user-provided syms function.
%        default value is 1e-8
%    .maxsteps    ... maximal number of integration steps. default is specified in the user-provided syms function.
%        default value is 1e4
%    .tstart    ... start of integration. for all timepoints before this, values will be set to initial value.
%        default value is 0
%    .sens_ind ... 1 dimensional vector of indexes for which sensitivities must be computed.
%        default value is 1:length(theta).
%    .x0 ... user-provided state initialisation. This should be a vector of dimension [#states, 1].
%        default is state initialisation based on the model definition.
%    .sx0 ... user-provided sensitivity initialisation. this should be a matrix of dimension [#states x #parameters].
%        default is sensitivity initialisation based on the derivative of the state initialisation.
%    .lmm    ... linear multistep method for forward problem.
%        1: Adams-Bashford
%        2: BDF (DEFAULT)
%    .iter    ... iteration method for linear multistep.
%        1: Functional
%        2: Newton (DEFAULT)
%    .linsol   ... linear solver module.
%        direct solvers:
%        1: Dense
%        2: Band (not implemented)
%        3: LAPACK Dense (not implemented)
%        4: LAPACK Band  (not implemented)
%        5: Diag (not implemented)
%        implicit krylov solvers:
%        6: SPGMR
%        7: SPBCG
%        8: SPTFQMR
%        sparse solvers:
%        9: KLU (DEFAULT)
%    .stldet   ... flag for stability limit detection. this should be turned on for stiff problems.
%        0: OFF
%        1: ON (DEFAULT)
%    .sensi   ... sensitivity order.
%        0: OFF (DEFAULT)
%        1: first
%        2: second
%    .sensi_meth   ... method for sensitivity analysis.
%        0: no sensitivity analysis
%        1 or 'forward': forward sensitivity analysis (DEFAULT)
%        2 or 'adjoint': adjoint sensitivity analysis 
%        3 or 'ss': defined but not documented 
%    .adjoint   ... flag for adjoint sensitivity analysis.
%        true: on 
%        false: off (DEFAULT)
%        NO LONGER USED: Replaced by sensi_meth
%    .ism   ... only available for sensi_meth == 1. Method for computation of forward sensitivities.
%        1: Simultaneous (DEFAULT)
%        2: Staggered
%        3: Staggered1
%    .Nd   ... only available for sensi_meth == 2. Number of Interpolation nodes for forward solution. 
%        default is 1000. 
%    .interpType   ... only available for sensi_meth == 2. Interpolation method for forward solution.
%        1: Hermite (DEFAULT for problems without discontinuities)
%        2: Polynomial (DEFAULT for problems with discontinuities)
%    .ordering   ... online state reordering.
%        0: AMD reordering (default)
%        1: COLAMD reordering
%        2: natural reordering
%    .newton_maxsteps   ... maximum newton steps
%        default value is 40
%        a value of 0 will disable the newton solver
%    .newton_maxlinsteps   ... maximum linear steps
%        default value is 100
%    .newton_preeq   ... preequilibration of system via newton solver
%        default value is false
%    .pscale   ... parameter scaling
%        []: (DEFAULT) use value specified in the model (fallback: 'lin')
%        0 or 'lin': linear
%        1 or 'log': natural log (base e)
%        2 or 'log10': log to the base 10
%
% Outputs:
% ========
% sol.status ... flag for status of integration. generally status<0 for failed integration
% sol.t ... vector at which the solution was computed
% sol.llh ... likelihood value
% sol.chi2 ... chi2 value
% sol.sllh ... gradient of likelihood
% sol.s2llh ... hessian or hessian-vector-product of likelihood
% sol.x ... time-resolved state vector
% sol.y ... time-resolved output vector
% sol.sx ... time-resolved state sensitivity vector
% sol.sy ... time-resolved output sensitivity vector
% sol.z ... event output
% sol.sz ... sensitivity of event output
function varargout = simulate_pd_fv(varargin)

% DO NOT CHANGE ANYTHING IN THIS FILE UNLESS YOU ARE VERY SURE ABOUT WHAT YOU ARE DOING
% MANUAL CHANGES TO THIS FILE CAN RESULT IN WRONG SOLUTIONS AND CRASHING OF MATLAB
if(nargin<2)
    error('Not enough input arguments.');
else
    tout=varargin{1};
    theta=varargin{2};
end
if(nargin>=3)
    kappa=varargin{3};
else
    kappa=[];
end

if(length(theta)<27)
    error('provided parameter vector is too short');
end


xscale = [];
if(nargin>=5)
    if(isa(varargin{5},'amioption'))
        options_ami = varargin{5};
    else
        options_ami = amioption(varargin{5});
    end
else
    options_ami = amioption();
end
if(isempty(options_ami.sens_ind))
    options_ami.sens_ind = 1:27;
end
if(options_ami.sensi>1)
    error('Second order sensitivities were requested but not computed');
end

if(isempty(options_ami.pscale))
    options_ami.pscale = 'log' ;
end
if(nargout>1)
    if(nargout>4)
        options_ami.sensi = 1;
        options_ami.sensi_meth = 'forward';
    else
        options_ami.sensi = 0;
    end
end
if(options_ami.sensi>0)
    if(options_ami.sensi_meth == 2)
        error('adjoint sensitivities are disabled as necessary routines were not compiled');
    end
end
nplist = length(options_ami.sens_ind); % MUST NOT CHANGE THIS VALUE
if(nplist == 0)
    options_ami.sensi = 0;
end
nxfull = 299;
plist = options_ami.sens_ind-1;
if(nargin>=4)
    if(isempty(varargin{4}))
        data=[];
    else
        if(isa(varargin{4},'amidata'))
             data=varargin{4};
        else
            data=amidata(varargin{4});
        end
        if(data.ne>0)
            options_ami.nmaxevent = data.ne;
        else
            data.ne = options_ami.nmaxevent;
        end
        if(isempty(kappa))
            kappa = data.condition;
        end
        if(isempty(tout))
            tout = data.t;
        end
    end
else
    data=[];
end
if(~all(tout==sort(tout)))
    error('Provided time vector is not monotonically increasing!');
end
if(max(options_ami.sens_ind)>27)
    error('Sensitivity index exceeds parameter dimension!')
end
if(length(kappa)<299)
    error('provided condition vector is too short');
end
init = struct();
if(~isempty(options_ami.x0))
    if(size(options_ami.x0,2)~=1)
        error('x0 field must be a column vector!');
    end
    if(size(options_ami.x0,1)~=nxfull)
        error('Number of rows in x0 field does not agree with number of states!');
    end
    init.x0 = options_ami.x0;
end
if(~isempty(options_ami.sx0))
    if(size(options_ami.sx0,2)~=nplist)
        error('Number of columns in sx0 field does not agree with number of model parameters!');
    end
    if(size(options_ami.sx0,1)~=nxfull)
        error('Number of rows in sx0 field does not agree with number of states!');
    end
end
sol = ami_pd_fv(tout,theta(1:27),kappa(1:299),options_ami,plist,xscale,init,data);
if(nargout>1)
    varargout{1} = sol.status;
    varargout{2} = sol.t;
    varargout{3} = sol.x;
    varargout{4} = sol.y;
    if(nargout>4)
        varargout{5} = sol.sx;
        varargout{6} = sol.sy;
    end
else
    varargout{1} = sol;
end
function chainRuleFactors = getChainRuleFactors(pscale, theta, sens_ind)
    if(length(pscale) == 1 && length(sens_ind) ~= length(pscale))
        chainRuleFactors = arrayfun(@(x, ip) getChainRuleFactor(x, theta(ip)), repmat(pscale, 1, length(sens_ind)), sens_ind);
    else
        chainRuleFactors = arrayfun(@(x, ip) getChainRuleFactor(x, theta(ip)), pscale, sens_ind);
    end
end

function chainRuleFactor = getChainRuleFactor(pscale, parameterValue)
    switch (pscale)
        case 1
            chainRuleFactor = exp(parameterValue);
        case 2
            chainRuleFactor = 10.^parameterValue*log(10);
        otherwise
            chainRuleFactor = 1.0;
    end
end

end
