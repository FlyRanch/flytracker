% GSSM_flyOcc  Template file for generalized state space model.

% The model is a 3-D Drosophila defined by 12 parameters, 6 rigid body 
% motion and 6 joint angles (3 for each wing).
%
% Process and observation noise are from Gaussian distribution.

%   This template file is used to completely describe a system in a generalized
%   state space format useable by the ReBEL inference and estimation system.
%   This file must be copied, renamed and adapted to your specific problem. The
%   interface to each function should NOT BE CHANGED however.
%
%   The following main and subfunctions must be defined:
%
%   1) [VARARGOUT] = MODEL_INTERFACE(FUNC, VARARGIN) :  This function is the
%        main gateway function which is used to initialise the generalized state
%        space model data structure. This is done by calling the 'init'
%        subfunction. The user can extend this function to indirectly call other
%        subfunctions within this file if needed.
%
%   2) MODEL = INIT(INIT_ARGS) : This function generates and initializes a
%        generalized state space model (gssm) data structure which summarizes
%        all relevant information about the system. 'model' is a Matlab structure
%        that must contain the following fields (consistentency can be checked
%        with the 'consistent' function.
%       - type : String which contains the model type. Use 'gssm' for generalized
%                state space model.
%       - tag  : ID string which contains instance specific identification info.
%                Default value=''
%       - ffun_type : State transition function type : 'lti' (linear time
%                     invariant), 'ltv' (linear time variant), 'nla' (nonlinear
%                     with additive noise) or 'nl' (pure nonlinear)
%       - hfun_type : State observation function type : 'lti', 'ltv', 'nla'
%                     or 'nl'.
%       - ffun : Function-handle to the state transition (state dynamics)
%                subfunction.
%       - hfun : Function-handle to the state observation subfunction.
%       - setparams : Function-handle to the setparams subfunction to update and
%                     possibly unpack the model parameters.
%
%       - prior : <<optional>> Function-handle to the state transition 'prior'
%                 function that calculates P(x(k)|x(k-1)). This must be defined
%                 if any of the particle filter family of estimators will be used
%                 on this model.
%       - likelihood : <<optional>> Function-handle to the observation likelihood
%                      function that calculates p(y(k)|x(k)) for a given
%                      realization of the state variable x and a particular
%                      observation instance y. This must be defined if any of the
%                      particle filter family of estimators will be used on this
%                      model.
%       - innovation : <<optional>> Function-handle to the innovation model
%                      function that calculates the difference between the output
%                      of the observation function (hfun) and the actual
%                      'real-world' measurement/observation of that signal. If
%                      this field is not defined, a generic innovation is used.
%       - linearize : <<optional>> Function-handle to the linearization
%                     subfunction. This is only needed if a linear Kalman filter
%                     (kf) or Extended Kalman Filter (ekf) will be used on this
%                     model. If no subfunction is defined, a default 'perturbation'
%                     based method of linearization will be used. This function
%                     does not need to be defined for the use of any of the
%                     Sigma-Point Kalman Filters (ukf, cdkf, srukf & srcdkf) or
%                     any of the Particle Filters (pf & sppf).
%
%       - statedim  : state dimension (this should be consistentent with the ffun
%                     and hfun subfunctions).
%       - obsdim    : observation dimension (this should be consistentent with
%                     the hfun subfunction).
%       - paramdim  : parameter dimension (number of free parameters in the
%                     system).
%       - U1dim     : dimension of exogenous input to ffun
%       - U2dim     : dimension of exogenous input to hfun
%       - pNoise    : process noise data structure  (this data structure is of
%                     type NoiseDS)
%       - oNoise    : observation noise data structure (this data structure is of
%                     type NoiseDS)
%       - params    : vector to hold all model parameters (must be of dimension
%                     [paramdim-by-1] )
%
%       - stateAngleCompIdxVec : <<optional>> Index vector idicating which (if
%                                any) of the state vector components are
%                                angular quantities (measured in radians) that
%                                has a discontinuety at +-Pi radians (this is
%                                needed by all SPKF based algorithms and derived
%                                hybrids)
%       - obsAngleCompIdxVec   : <<optional>> Index vector idicating which (if
%                                any) of the observation vector components are
%                                angular quantities (measured in radians) that
%                                has a discontinuety at +-Pi radians (this is
%                                needed by all SPKF based algorithms and
%                                derived hybrids)
%
%   3) MODEL = SETPARAMS(MODEL, PARAMS, IDXVECTOR) : This function unpacks a
%      column vector containing system parameters into specific forms needed by
%      FFUN, HFUN and possibly defined sub-functional objects. Both the
%      vectorized (packed) form of the parameters 'PARAMS' as well as the
%      unpacked forms are stored within the model data structure. 'IDXVECTOR' is
%      an optional argument which indicates which parameters should be updated.
%      This can be used to only modify a subset of the total system parameters.
%      'PARAMS' and 'IDXVECTOR' must have the same length.
%      Example :  model=setparams(model, [1 1 2 1 3]', [1 3 6:8])
%      << THIS SUBFUNCTION IS REQUIRED >>
%
%   4) NEW_STATE = FFUN(MODEL, STATE, V, U1) : State transition function which
%      takes as input the current state of the system 'STATE', a process noise
%      vector 'V', an exogenous control input 'U1', and a gssm data structure
%      'MODEL', and calculates the system state at the next discrete time instant,
%      'NEW_STATE'. This function implements the system dynamics.
%      << THIS SUBFUNCTION IS REQUIRED >>
%
%   5) OBSERV = HFUN(MODEL, STATE, N, U2) : State observation function which
%      takes as input the current state of the system 'STATE', an observation
%      noise vector 'N', an exogenous control input 'U2' and a gssm data
%      structure 'MODEL', and calculates the current observation vector of the
%      system, 'OBSERV'.
%      << THIS SUBFUNCTION IS REQUIRED >>
%
%   6) TRAN_PRIOR = PRIOR(MODEL, NEXT_STATE, STATE, U1, PNOISEDS) : Calculates
%      the transition prior p(next_state|state) = p(state(k)|state(k-1)) =
%      p(x(k)|x(k-1)) given a gssm data structure 'MODEL', realizations of the
%      system state at time k and k-1, 'NEXT_STATE' and 'STATE' and the
%      exogeneous inputs to the process model, U1. The process noise data
%      structure 'PNOISEDS' specifies which noise model should be used to
%      calculate the likelihood. If this is ommitted, the default model defined
%      process noise data structure 'model.pNoise' is used.
%      << THIS SUBFUNCTION IS OPTIONAL : Only required by particle filter
%         family of estimator >>
%
%   7) LLH = LIKELIHOOD(MODEL, OBS, STATE, U2, ONOISEDS) : Calculates the
%      likelihood of a 'real world' observation 'OBS' for a given realization
%      or instance of the state variable STATE. i.e. Calculates the value of
%      P(OBS|STATE). The measurement noise data structure 'ONOISEDS' specifies
%      which noise model should be used to calculate the likelihood. If this is
%      ommitted, the default model defined observation noise data structure
%      'model.pNoise' is used. 'U2' is the (optional) exogeneous input to the
%      state observation function 'hfun'.
%      << THIS SUBFUNCTION IS OPTIONAL : Only required by particle filter
%         family of estimator >>
%
%   8) INNOV = INNOVATION(MODEL, OBS, OBSERV) : Calculates the innovation signal
%      (difference) between the output of HFUN, i.e. OBSERV=HFUN(STATE) (the
%      predicted system observation) and an actual 'real world' observation OBS.
%      This function might be as simple as INNOV = OBS - OBSERV, which is the
%      default case, but can also be more complex for complex measurement
%      processes where for example multiple (possibly false) observations can be
%      observed for a given hidden ground truth.
%      << THIS SUBFUNCTION IS OPTIONAL : Only redefine if the default does not
%         reflect the true measurement process >>
%
%   9) OUT = LINEARIZE(MODEL, STATE, V, N, U1, U2, TERM, IDXVECTOR) generates a
%      linearized model of the nonlinear system described by the gssm data
%      structure MODEL at the current operating point, STATE, exogenous inputs U1
%      and U2. The linearized model is of the form:
%
%           state(k) = A*state(k-1) + B*u1(k-1) + G*v(k-1)
%               y(k) = C*state(k)   + D*u2(k)   + H*n(k)
%
%      for an arbitrary model defined by this GSSM file. The string TERM
%      specifies which of the model terms are returned, i.e.
%
%       A = linearize(model, state, v, n, u1, u2, 'A') or
%       H = linearize(model, state, v, n, u1, u2, 'H') etc.
%
%      TERM can be one of the following, 'A','B','C','D','G','H','JFW','JHW' ,
%      where 'JFW' and 'JHW' are the partial derivatives of FFUN and HFUN with
%      respect to the system parameters.
%
%      IDXVECTOR is an optional argument indicating which subset of the
%      independent vector should be used to calculate any specific derivative.
%      This will result in a Jacobian matrix with a reduced number of columns,
%      corresponding with the subvector as defined by 'IDXVECTOR'. The default
%      (when this argument is ommitted) is to use the full vector.
%
%      << THIS SUBFUNCTION IS OPTIONAL : Only required for Kalman and Extended
%      Kalman filters >>
%
%     See also
%     CONSIST, GENINFDS, GENNOISEDS
%

%   Copyright  (c) Rudolph van der Merwe (2002)
%
%   This file is part of the ReBEL Toolkit. The ReBEL Toolkit is available free for
%   academic use only (see included license file) and can be obtained by contacting
%   rvdmerwe@ece.ogi.edu.  Businesses wishing to obtain a copy of the software should
%   contact ericwan@ece.ogi.edu for commercial licensing information.
%
%   See LICENSE.TXT (which should be part of the main toolkit distribution) for more
%   detail.


%===============================================================================================

function [varargout] = model_interface(func, varargin)
global PAR
  switch func

    %--- Initialize GSSM data structure --------------------------------------------------------
    case 'init'
      model = init(varargin);
      error(consistent(model,'gssm'));              % check consistentency of initialized model
      varargout{1} = model;


      %--------------------------------------------------------------------------------------------
    otherwise

      error(['Function ''' func ''' not supported.']);

  end


%===============================================================================================
function model = init(init_args)

global PAR

  model.type = 'gssm';                       % object type = generalized state space model

  model.tag  = 'GSSM_fly12';                           % ID tag

  model.ffun_type = 'nla';                   % FFUN type     : 'lti'  - linear time invariant
                                             %                 'ltv'  - linear time variant
                                             %                 'nla'  - nonlinear with additive noise
                                             %                 'nl'   - nonlinear

  model.hfun_type = 'nla';                   % HFUN type     : 'lti', 'ltv', 'nla' or 'nl'

  model.setparams = @setparams;              % function handle to SETPARAMS
  model.ffun      = @ffun;                   % function handle to FFUN
  model.hfun      = @hfun;                   % function handle to HFUN
  model.cfun      = @cfun; 
  % model.prior = @prior;                    % function handle to PRIOR        (uncomment if 'prior' subfunction is defined)
  % model.likelihood = @likelihood;          % function handle to LIKELIHOOD   (uncomment if 'likelihood' subfunction is defined)
  % model.innovation = @innovation;          % function handle to INNOVATION   (uncomment if 'innovation' subfunction is defined)
  % model.linearize = @linearize;              % function handle to LINEARIZE    (uncomment if 'linearize' subfunction is defined)

  % model.stateAngleCompIdxVec :             % <<optional>> Index vector idicating which (if
                                             % any) of the state vector components are
                                             % angular quantities (measured in radians) that
                                             % has a discontinuety at +-Pi radians (this is
                                             % needed by all SPKF based algorithms and derived
                                             % hybrids)
  % model.obsAngleCompIdxVec                 % <<optional>> Same as model.stateAngleCompIdxVec but
                                             % used for the observation vector.

  model.statedim   = PAR.numfly*PAR.statedim;      % state dimension
  % L = number of sample points along the length of fly
  % D = dimension of each observation (1 for normal projection, 2
  %                                    for euclidean distance)
  L = PAR.L1;
  D1 = 1;
  D2 = 2;
  Npix = 0;
  
  % The observation dimension per fly is:
  % L*2sides*D + 2HeadTail*D
  %model.obsdimPER  = L*2*D2 + 1*D2; %for modcurve8newP  % observation dimension per model 
  model.obsdimPER  = 2*L*D1 + 1*D1; %for modcurve8newP1
%   model.obsdimPER  = 2*L*D1 + 1*D1 + 1*D2; %for modcurve8newP2
  model.obsdim     = model.obsdimPER*PAR.numfly;   % observation dimension
  
  model.Pixobsdim = Npix;
  model.paramdim  = PAR.paramdim*PAR.numfly;      % total parameter dimension
                                                  % (wormradius,length)
  model.U1dim     = 0;                            % exogenous control input 1 dimension
  model.U2dim     = 0;                            % exogenous control input 2 dimension
  model.Vdim      = PAR.pNoisedim*PAR.numfly;     % process noise dimension
  % we want to add noise to each component of the observation
  % feature points.
  model.Ndim       = model.obsdimPER*PAR.numfly;   % observation noise dimension

  %-----------------------------------------
  %-- Setup process noise source
  Arg.type = 'gaussian';                     % noise source type : Gaussian
  Arg.cov_type = 'full';                     % Gaussian noise source cov_type (full covariance)
  Arg.tag = 'GSSM process noise source';     % Arbitrary ID tag (optional)
  Arg.dim = model.Vdim;                      % noise dimension
  Arg.mu = zeros(model.Vdim,1);              % noise mean
  
  
  %This is the uncertainty (2*sigma) in degrees of joint angles for the Fly
  %wing
  twoSigmaDeg = [10];
  WingVar = (twoSigmaDeg.*(pi/180).*.5).^2;
  
  %Variance for the wing joint locations
  JointVar = [0.001 0.0004 0.0003];
  
  %Variance for the body linear acceleration
  LinVar = [3 2 4]*1e11;
  %Variance for the body angular acceleration
  AngVar = [2 23 3]*1e11;
  
  %Variance for the body rollangular acceleration
  RollVar = 4e9;
  
  %fly Process Noise
  Qbody = [0.001.*ones(1,3) 0.00001.*ones(1,3)];
  Q1 = [Qbody repmat(WingVar,1,6)];
  %Q1 = [Qbody repmat(WingVar,1,6) repmat(JointVar,1,2)];
  %Q1 = [repmat(WingVar,1,6) LinVar AngVar RollVar];
  
  Q = diag(repmat(Q1,1,PAR.numfly));
  Arg.cov  = Q;            % noise covariance
  model.pNoise = gennoiseds(Arg);            % generate noise source

  %-----------------------------------------
  %-- Setup observation noise source
  Arg.type = 'gaussian';                     % noise source type : Gaussian
  Arg.cov_type = 'diag';                      % Gaussian noise source cov_type (full covariance)
  Arg.tag = 'GSSM observation noise source'; % Arbitrary ID tag (optional)
  Arg.dim = model.Ndim;                      % noise dimension
  Arg.mu = zeros(model.Ndim,1);              % noise mean
  
  % This corresponds to a variance of 1 pixel for each boundary
  % point that is observed
%   P = [3 3 .25 .25 .25 0.5*ones(1,model.obsdimPER - 5)]; %for modcurve8newP2
  P = [.25 .25 .25 0.5*ones(1,model.obsdimPER - 3)]; %for modcurve8newP1 
  %P3 = [1 1 1 4e-2*ones(1,3) ones(1,model.obsdimPER-6)]; 
%   Arg.cov = P3;  % noise covariance
%   Arg.cov = blkdiag(Arg.cov,diag(repmat(P,1,PAR.numfly)));  % noise covariance
  Arg.cov = diag(repmat(P,1,PAR.numfly));
%   Arg.cov = diag(repmat(P3,1,PAR.numfly));  % noise covariance
    
  model.oNoise = gennoiseds(Arg);            % generate noise source
    
  model.params     = zeros(model.paramdim,1); %  setup parameter vector buffer
  % params is a structure with all the fly geometric data in it. 
  model = setparams(model,PAR.params);   % initialize model parameters and unpack if needed
  
%-- The subsection of the init function below these comments should be used to define any other data structures, objects,
%-- functions, etc.which is needed by the internal implementation of the FFUN, HFUN, LINEARIZE, ETC. functions. These data
%-- structures should be saved within the GSSM 'model' data structure. The user can embed any other structures such as Netlab
%-- neural networks, etc. in this section.


%===============================================================================================
function model = setparams(model, params, idxVector)

% Function to unpack a column vector containing system parameters into specific forms
% needed by FFUN, HFUN and possibly defined sub-functional objects. Both the vectorized (packed)
% form of the parameters as well as the unpacked forms are stored within the model data structure.
% INDEX_VECTOR is an optional argument which indicates which parameters should be updated. This can
% be used to only modify a subset of the total system parameters.
%
%   Exmaple: model=setparams(model, [1 1 2 1 3]', [1 3 6:8]);

  switch nargin

  case 2
    %------------  Set all system parameters ---------------------------------------------------
    
    model.params = params;
    
    %-- Add unpack code here if needed -----
    %---------------------------------------
    
    
   case 3
    %------------  Set a subset of system parameters -------------------------------------------
    
    model.params(idxVector) = params;
    
    %-- Add unpack code if needed here.
    %---------------------------------------
    
   otherwise
    error('[ setparams ] Incorrect number of input arguments.');
    
  end
  
  
%===============================================================================================
function [new_state,model] = ffun(model, state, V, U1)

% FFUN  State transition function (system dynamics).
%
%   Generates the next state of the system NEW_STATE given
%   the current STATE, exogenous input U1 and process noise term V. If STATE, U1 and V are matrices
%   then FFUN is calculated for each column vector of these matrices, resulting in an equal number
%   of columns in NEW_STATE. MODEL is a GSSM derived data structure describing the system
global PAR

%---Initialize
new_state = zeros(size(state));

%Right now, just do a random walk in parameters
%keyboard
if ~isempty(V)
    % Additive noise for Shape parameters and Translation
    for k = 1:PAR.numfly
        new_state = state + V;
    end
end

    


%===============================================================================================
function  tranprior = prior(model, nextstate, state, U1, pNoiseDS)

% PRIOR  Transition prior function
%
%   Calculates P(nextstate|state). If you plan to run a particle filter on this mode, you should
%   define this.
%
%   INPUT
%         model          GSSM data structure
%         nextstate      state at time k
%         state          state at time k-1
%         U1             exogeneous input to FFUN at time k-1
%         pNoiseDS       (optional) process noise NoiseDS data structure to use for evaluation of
%                        transition prior. If this is ommitted, model.pNoise, is used.
%   OUTPUT
%         tranprior      p(x(k)|x(k-1))
%
%-- This function must be defined by the user!

%===============================================================================================
function [observ,model] = hfun(model, state, N, U2)

% HFUN  State observation function.
%
%   OBSERV = HFUN(MODEL, STATE, N, U2) generates the current possibly nonlinear observation of the
%   system state, OBSERV, given the current STATE, exogenous input U and observation noise term V.
%   If STATE, U2 and N are matrices then HFUN is calculated for each column vector of these matrices,
%   resulting in an equal number of columns in OBSERV. MODEL is a GSSM derived data structure describing
%   the system.
global PAR

%skew matrix to calculate cross-product
skew = @(p) [0 -p(3) p(2);p(3) 0 -p(1);-p(2) p(1) 0];

mdlLen = PAR.statedim;

% This observation function outputs the cross product of the 3D model
% points with the vector representing the direction of their projection
% ray
%--- Initialize

observ = zeros(model.obsdim,size(state,2));

for j = 1:size(state,2)
    %Iterate over the state vector for each Fly.

    a = [];
    for k = 1:PAR.numfly
         pp = state((k-1)*mdlLen+1:k*mdlLen,j);
         
         %Evaluate model at current state
         [x,y,z] = flymod(pp,PAR.params,PAR);
         %[x,y,z] = flymod_BodyOnly(pp,PAR.params,PAR);

         for i = 1:length(x)
             PAR.modsample(i) = size(x{i},1);
             % reshape surface matrix structure to N x 3 matrix of all points
             % for ith part
             pts{i} = [reshape(x{i},[],1) reshape(y{i},[],1) reshape(z{i},[],1)];
         end
         
         % Now iterate over each camera
         for cam = 1:PAR.numcam
             %Get the points that are corresponding to the 2D boundary
             pts3D = [];
             %iterate over each part
             for i = 1:length(pts)
                 index = model.IdxTo3DPts{cam,k}{i};
                 temp = pts{i}(index,:);
                 pts3D = [pts3D ; temp];
             end
             
             % Get rid of the points that are occluded
             occidx = model.occluded_idx{cam,k};
             pts3D(occidx,:) = [];
             
             
             
             %Now I want to take the cross product of these points with the
             %projection rays
             
             %line directions: First 3 components of ray's Plüker
             %coordinates
             Direction = model.DataRays{cam}(:,1:3);
             
             %Let's create a giant sparse skew-symmetric matrix to perform
             %the cross product of all these vectors with 1 matrix
             %operation!
             nn = size(pts3D,1);
             % The code below works but the addition for BigSkewMatrix is
             % too slow!
             % Instead create sparse matrix directly
%              upperZpart = sparse(1:3:3*nn,2:3:3*nn,-pts3D(:,3),3*nn,3*nn);
%                          
%              upperYpart = sparse(1:3:3*nn,3:3:3*nn,pts3D(:,2),3*nn,3*nn);
%              
%              upperXpart = sparse(2:3:3*nn,3:3:3*nn,-pts3D(:,1),3*nn,3*nn);
%              
%              % These parts create the upper triangle of skew matrix
%              % Transpose them and change the sign to create the lower
%              % triangle
%              BigSkewMatrix = upperZpart + upperYpart + upperXpart -... 
%                  (upperZpart' + upperYpart' + upperXpart');  
              rows = [1:3:3*nn 2:3:3*nn 1:3:3*nn 3:3:3*nn 2:3:3*nn 3:3:3*nn];
              cols = [2:3:3*nn 1:3:3*nn 3:3:3*nn 1:3:3*nn 3:3:3*nn 2:3:3*nn];
              ppts = [-pts3D(:,3) ; pts3D(:,3) ; pts3D(:,2) ; -pts3D(:,2) ; -pts3D(:,1) ; pts3D(:,1)];
              BigSkewMatrix = sparse(rows,cols,ppts,3*nn,3*nn);
             
             
             crossprod = BigSkewMatrix*reshape(Direction',[],1);
             %crossprod = reshape(crossprod,3,[])';
             
%              crossprod = zeros(size(pts3D));
%              for i = 1:size(pts3D,1)
%                  crossprod(i,:) = ( ([0 -pts3D(i,3) pts3D(i,2)
%                                       pts3D(i,3) 0 -pts3D(i,1)
%                                       -pts3D(i,2) pts3D(i,1) 0])*Direction(i,:)' )';
%              end
             
             %a = [a ; reshape(crossprod',[],1)];
             
             % Calculate the scaling to make it be the 2D image projection error
             % See B. Rosenhahn, et. al., "Three-Dimensional Shape Knowledge for
             % Joint Image Segmentation and Pose Tracking" IJCV, 2007.
             % eqn. 3.20
%              Moments = model.DataRays{cam}(:,4:6);
%              obs = reshape(Moments',[],1);
%              
%              inov = crossprod - obs;
%              inov_mat = reshape(inov,3,[]);
%              eta_denom = sqrt( sum(inov_mat.^2,1))';
%              eta_scale = model.eta_Num{cam,k} ./ eta_denom;
%              % I need to take the sqrt to match his equation.  Also, repeat scale
%              % 3 times for each point for the 3D coordinates
%              eta_scale = reshape( repmat(eta_scale',3,1),[],1);
%              eta_scale = sqrt(eta_scale);
%              
%              %Store the scaling value
%              model.eta_scale{cam,k} = eta_scale;
             
             a = [a ; crossprod];
             %keyboard
         end
    end
    observ(:,j) = a;
end

if ~isempty(N),
    observ = observ + N;
end


% %===============================================================================================
% function [new_state,model] = cfun(model, state)
% 
% % CFUN  constraint function .
% %
% %   Generates the next state of the system NEW_STATE given
% %   the current STATE, exogenous input U1 and process noise term V. If STATE, U1 and V are matrices
% %   then FFUN is calculated for each column vector of these matrices, resulting in an equal number
% %   of columns in NEW_STATE. MODEL is a GSSM derived data structure describing the system
% %
% % Constrain the joint location displacements to be less than 'const'
% global PAR
% 
% %These are the indices of the state that have a constraint associated with
% %them, and the constraint value
% cidx = {13:15,16:18};
% const = {0.09,0.09};
% 
% %---Initialize
% new_state = zeros(size(state));
% 
% mdlLen = PAR.statedim;
% 
% for k = 1:PAR.numfly
%     xh = state((k-1)*mdlLen+1:k*mdlLen,:);
%     %Apply Constraints
%     for j = 1:length(cidx)
%         %Calculate norm
%         NRM = sqrt(sum(xh(cidx{j},:).^2,1));
%         %Rescale the values that have norm greater that const
%         if any(NRM > const{j})
%             xh(cidx{j},NRM > const{j}) = const{j}.*xh(cidx{j},NRM > const{j})...
%                 ./repmat(NRM(NRM > const{j}),length(cidx{j}),1);
%         end
%     end
% 
%     %Predict the new state parameters
%     new_state((k-1)*mdlLen+1:k*mdlLen,:) = xh;
% end

%===============================================================================================
function [new_state,model] = cfun(model, state)

% CFUN  constraint function .
%
%   Generates the next state of the system NEW_STATE given
%   the current STATE, exogenous input U1 and process noise term V. If STATE, U1 and V are matrices
%   then FFUN is calculated for each column vector of these matrices, resulting in an equal number
%   of columns in NEW_STATE. MODEL is a GSSM derived data structure describing the system
%
% Constrain the roll angle so that it is symmetric between the wing joints
% in the transverse plane

global PAR

%---Initialize
new_state = zeros(size(state));

mdlLen = PAR.statedim;

for i = 1:size(state,2)
    for k = 1:PAR.numfly
        xh = state((k-1)*mdlLen+1:k*mdlLen,i);
        %Apply Constraints
        new_state((k-1)*mdlLen+1:k*mdlLen,i) = fixRollTPlane(xh,PAR);
    end
end

%==========================================================================


%===============================================================================================
function llh = likelihood(model, obs, state, U2, oNoiseDS)

% LIKELIHOOD  Observation likelihood function
%
% Function-handle to the observation likelihood function that calculates p(y|x) for a
% given realization of the state variable 'state' and a particular observation instance 'obs'.
%
%   i.e. Calculates the value of P(OBS|STATE) = P(y|x)
%
%   INPUT
%         model          GSSM data structure
%         obs            observation at time k
%         state          state at time k
%         U2             exogeneous input to HFUN at time k
%         oNoiseDS       (optional) measurement noise NoiseDS data structure to use for evaluation of
%                        transition prior. If this is ommitted, model.oNoise, is used.
%   OUTPUT
%         llh            p(y(k)|x(k))
%
%-- This function must be defined by the user!


%===============================================================================================
function INNOV = innovation(model, obs, observ)

% INNOVATION  Innovation model
%
%   INNOV = INNOVATION(MODEL, STATE, OBS, OBSERV) : Calculates the innovation signal (difference) between the
%   output of HFUN, i.e. OBSERV (the predicted system observation) and an actual 'real world' observation OBS.
%   This function might be as simple as INNOV = OBS - OBSERV, which is the default case, but can also be more
%   complex for complex measurement processes where for example multiple (possibly false) observations can be
%   observed for a given hidden ground truth.

%-- This function must be redefined by the user if the specific real world observation process dictates it

%-- Acquire the normal vectors from the data.


%===============================================================================================
function out = linearize(model, state, V, N, U1, U2, term, idxVector)

% LINEARIZE
%
%   OUT = LINEARIZE(MODEL, STATE, V, N, U1, U2, TERM, IDXVECTOR) returns a linearized model of the
%   form
%           state(k) = A*state(k-1) + B*u1(k-1) + G*v(k-1)
%               y(k) = C*state(k)   + D*u2(k)   + H*n(k)
%
%   for an arbitrary model defined by this GSSM file. The string TERM specifies which of the
%   model terms are returned, i.e.
%
%   A = linearize(model, state, v, n, u1, u2, 'A') or
%   O = linearize(model, state, v, n, u1, u2, 'H') etc.
%
%   TERM can be one of the following, 'A','B','C','D','G','H','JFW','JHW' , where 'JFW' and 'JHW'
%   are the partial derivatives of FFUN and HFUN with respect to the system parameters.
%
%   INDEX_VECTOR is an optional argument indicating which subset of the independent vector should be used to calculate
%   any specific derivative. This will result in a Jacobian matrix with a reduced number of columns,
%   corresponding with the subvector as defined by index_vector. The default (when this argument is ommitted)
%   is to use the full vector.
%
%   Generic perturbation based linearization subunits are provided. These can (and should) be replaced
%   by user defined analytical derivative code if available. If no linearization function is available
%   or is not needed, a call to this function should return an error message.
%

  nia = nargin;                      % number of input arguments
  if (nia < 7)
    error('[ linearize ] Not enough input arguments! ');
  end

  epsilon = 1e-8;                    % perturbation step size

  switch (term)

    case 'A'
      %%%========================================================
      %%%             Calculate A = dffun/dstate
      %%%========================================================
      if (nia==7), index_vector=[1:model.statedim]; end
      liv = length(index_vector);
      A  = zeros(model.statedim, liv);
      %%%---------- replace this section if needed --------------
      f1 = feval(model.ffun,model,state,V,U1);
      for j=1:liv,
        s = state;
        k = index_vector(j);
        s(k) = s(k) + epsilon;
        f2 = feval(model.ffun,model,s,V,U1);
        A(:,j) = (f2-f1)/epsilon;
      end
      %%%--------------------------------------------------------
      out = A;


    case 'B'
      %%%========================================================
      %%%             Calculate B = dffun/dU1
      %%%========================================================
      if (nia==7), index_vector=[1:model.U1dim]; end
      liv = length(index_vector);
      B = zeros(model.statedim, liv);
      %%%---------- replace this section if needed --------------
      f1 = feval(model.ffun,model,state,V,U1);
      for j=1:liv,
        Utemp = U1;
        k = index_vector(j);
        Utemp(k) = Utemp(k) + epsilon;
        f2 = feval(model.ffun,model,state,V,Utemp);
        B(:,j) = (f2-f1)/epsilon;
      end
      %%%--------------------------------------------------------
      out = B;


    case 'C'
      %%%========================================================
      %%%             Calculate C = dhfun/dx
      %%%========================================================
      if (nia==7), index_vector=[1:model.statedim]; end
      liv = length(index_vector);
      C = zeros(model.obsdim, liv);
      %%%---------- replace this section if needed --------------
      f3 = feval(model.hfun,model,state,N,U2);
      for j=1:liv,
        s = state;
        k = index_vector(j);
        s(k) = s(k) + epsilon;
        f4 = feval(model.hfun,model,s,N,U2);
        C(:,j) = (f4-f3)/epsilon;
      end
      %%%--------------------------------------------------------
      out = C;


    case 'D'
      %%%========================================================
      %%%             Calculate D = dhfun/dU2
      %%%========================================================
      if (nia==7), index_vector=[1:model.U2dim]; end
      liv = length(index_vector);
      D = zeros(model.obsdim, liv);
      %%%---------- replace this section if needed --------------
      f3 = feval(model.hfun,model,state,N,U2);
      for j=1:liv,
        Utemp = U2;
        k = index_vector(j);
        Utemp(k) = Utemp(k) + epsilon;
        f4 = feval(model.hfun,model,state,N,Utemp);
        D(:,j) = (f4-f3)/epsilon;
      end
      %%%--------------------------------------------------------
      out = D;


    case 'G'
      %%%========================================================
      %%%             Calculate G = dffun/dv
      %%%========================================================
      if (nia==7), index_vector=[1:model.Vdim]; end
      liv = length(index_vector);
      G = zeros(model.statedim, liv);
      %%%---------- replace this section if needed --------------
      f1 = feval(model.ffun,model,state,V,U1);
      for j=1:liv,
        Vtemp = V;
        k = index_vector(j);
        Vtemp(k) = Vtemp(k) + epsilon;
        f5 = feval(model.ffun,model,state,Vtemp,U1);
        G(:,j) = (f5-f1)/epsilon;
      end
      %%%--------------------------------------------------------
      out = G;


    case 'H'
      %%%========================================================
      %%%             Calculate H = dhfun/dn
      %%%========================================================
      if (nia==7), index_vector=[1:model.Ndim]; end
      liv = length(index_vector);
      H = zeros(model.obsdim, liv);
      %%%---------- replace this section if needed --------------
      f3 = feval(model.hfun,model,state,N,U2);
      for j=1:liv,
        Ntemp = N;
        k = index_vector(j);
        Ntemp(k) = Ntemp(k) + epsilon;
        f6 = feval(model.hfun,model,state,Ntemp,U2);
        H(:,j) = (f6-f3)/epsilon;
      end
      %%%--------------------------------------------------------
      out = H;


    case 'JFW'
      %%%========================================================
      %%%             Calculate  = dffun/dparameters
      %%%========================================================
      if (nia==7), index_vector=[1:model.paramdim]; end
      liv = length(index_vector);
      JFW = zeros(model.statedim, liv);
      %%%---------- replace this section if needed --------------
      f1 = feval(model.ffun,model,state,V,U1);
      old_params = model.params;                         % save current model parameters
      for j=1:liv,
        params = old_params;
        k = index_vector(j);
        params(k) = params(k) + epsilon;
        model = setparams(model,params);
        f7 = feval(model.ffun,model,state,V,U1);
        JFW(:,j) = (f7-f1)/epsilon;
      end
      %%%--------------------------------------------------------
      out = JFW;


    case 'JHW'
      %%%========================================================
      %%%             Calculate  = dhfun/dparameters
      %%%========================================================
      if (nia==7), index_vector=[1:model.paramdim]; end
      liv = length(index_vector);
      JHW = zeros(model.obsdim, liv);
      %%%---------- replace this section if needed --------------
      f3 = feval(model.hfun,model,state,N,U2);
      old_params = model.params;                         % save current model parameters
      for j=1:liv,
        params = old_params;
        k = index_vector(j);
        params(k) = params(k) + epsilon;
        model = setparams(model,params);
        f8 = feval(model.hfun,model,state,N,U2);
        JHW(:,j) = (f8-f3)/epsilon;
      end
      %%%--------------------------------------------------------
      out = JHW;

    otherwise
      error('[ linearize ] Invalid linearization term requested!');

  end

  %--------------------------------------------------------------------------------------
