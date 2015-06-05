function results = analyzePRF_cmrr(userid,subject,stimulus,data,tr,options)

% function results = analyzePRF_chpc(userid,subject,stimulus,data,tr,options)
%
% the documentation of this function is identical to that for analyzePRF.m,
% except that we take two additional inputs:
%
% <userid> is a string with the user ID (this is used for submitting jobs)
% <subject> is a positive integer
%
% history:
% - 2015/02/07 - inherit from analyzePRF.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REPORT

fprintf('*** analyzePRF_chpc: started at %s. ***\n',datestr(now));
stime = clock;  % start time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTERNAL CONSTANTS

% define
scratchdir = sprintf('/tmp/%s/analyzePRF',userid);
mkdirquiet(scratchdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP AND PREPARATION

% massage cell inputs
if ~iscell(stimulus)
  stimulus = {stimulus};
end
if ~iscell(data)
  data = {data};
end

% calc
is3d = size(data{1},4) > 1;
if is3d
  dimdata = 3;
  dimtime = 4;
  xyzsize = sizefull(data{1},3);
else
  dimdata = 1;
  dimtime = 2;
  xyzsize = size(data{1},1);
end
numvxs = prod(xyzsize);

% calc
res = sizefull(stimulus{1},2);
resmx = max(res);
numruns = length(data);

% deal with inputs
if ~exist('options','var') || isempty(options)
  options = struct();
end
if ~isfield(options,'vxs') || isempty(options.vxs)
  options.vxs = 1:numvxs;
end
if ~isfield(options,'wantglmdenoise') || isempty(options.wantglmdenoise)
  options.wantglmdenoise = 0;
end
if ~isfield(options,'hrf') || isempty(options.hrf)
  options.hrf = [];
end
if ~isfield(options,'maxpolydeg') || isempty(options.maxpolydeg)
  options.maxpolydeg = [];
end
if ~isfield(options,'seedmode') || isempty(options.seedmode)
  options.seedmode = [0 1 2];
end
if ~isfield(options,'xvalmode') || isempty(options.xvalmode)
  options.xvalmode = 0;
end
if ~isfield(options,'numperjob') || isempty(options.numperjob)
  options.numperjob = [];
end
if ~isfield(options,'maxiter') || isempty(options.maxiter)
  options.maxiter = 500;
end
if ~isfield(options,'display') || isempty(options.display)
  options.display = 'iter';
end
if ~isfield(options,'typicalgain') || isempty(options.typicalgain)
  options.typicalgain = 10;
end

% massage
wantquick = isequal(options.seedmode,-2);
options.seedmode = union(options.seedmode(:),[]);

% massage more
if wantquick
  opt.xvalmode = 0;
  opt.vxs = 1:numvxs;
  opt.numperjob = [];
end

% calc
usecluster = ~isempty(options.numperjob);

% prepare stimuli
for p=1:length(stimulus)
  stimulus{p} = squish(stimulus{p},2)';  % frames x pixels
  stimulus{p} = [stimulus{p} p*ones(size(stimulus{p},1),1)];  % add a dummy column to indicate run breaks
  stimulus{p} = single(stimulus{p});  % make single to save memory
end

% deal with data badness (set bad voxels to be always all 0)
bad = cellfun(@(x) any(~isfinite(x),dimtime) | all(x==0,dimtime),data,'UniformOutput',0);  % if non-finite or all 0
bad = any(cat(dimtime,bad{:}),dimtime);  % badness in ANY run
for p=1:numruns
  data{p}(repmat(bad,[ones(1,dimdata) size(data{p},dimtime)])) = 0;
end

% calc mean volume
meanvol = mean(catcell(dimtime,data),dimtime);

% what HRF should we use?
if isempty(options.hrf)
  options.hrf = getcanonicalhrf(tr,tr)';
end
numinhrf = length(options.hrf);

% what polynomials should we use?
if isempty(options.maxpolydeg)
  options.maxpolydeg = cellfun(@(x) round(size(x,dimtime)*tr/60/2),data);
end
if isscalar(options.maxpolydeg)
  options.maxpolydeg = repmat(options.maxpolydeg,[1 numruns]);
end
fprintf('using the following maximum polynomial degrees: %s\n',mat2str(options.maxpolydeg));

% initialize cluster stuff
if usecluster
  filestodelete = {};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE OUT NOISE REGRESSORS

if isequal(options.wantglmdenoise,1)
  noisereg = analyzePRFcomputeGLMdenoiseregressors(stimulus,data,tr);
elseif isequal(options.wantglmdenoise,0)
  noisereg = [];
else
  noisereg = options.wantglmdenoise;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE MODEL

% pre-compute some cache
[d,xx,yy] = makegaussian2d(resmx,2,2,2,2);

% define the model (parameters are R C S G N)
modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),options.hrf,dd(:,prod(res)+1));
model = {{[] [1-res(1)+1 1-res(2)+1 0    0   NaN;
              2*res(1)-1 2*res(2)-1 Inf  Inf Inf] modelfun} ...
         {@(ss)ss [1-res(1)+1 1-res(2)+1 0    0   0;
                   2*res(1)-1 2*res(2)-1 Inf  Inf Inf] @(ss)modelfun}};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE SEEDS

% init
seeds = [];

% generic large seed
if ismember(0,options.seedmode)
  seeds = [seeds;
           (1+res(1))/2 (1+res(2))/2 resmx/4*sqrt(0.5) options.typicalgain 0.5];
end

% generic small seed
if ismember(1,options.seedmode)
  seeds = [seeds;
           (1+res(1))/2 (1+res(2))/2 resmx/4*sqrt(0.5)/10 options.typicalgain 0.5];
end

% super-grid seed
if any(ismember([2 -2],options.seedmode))
  [supergridseeds,rvalues] = analyzePRFcomputesupergridseeds(res,stimulus,data,modelfun, ...
                                                   options.maxpolydeg,dimdata,dimtime, ...
                                                   options.typicalgain,noisereg);
end

% make a function that individualizes the seeds
if exist('supergridseeds','var')
  seedfun = @(vx) [[seeds];
                   [subscript(squish(supergridseeds,dimdata),{vx ':'})]];
else
  seedfun = @(vx) [seeds];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERFORM OPTIMIZATION

% if this is true, we can bypass all of the optimization stuff!
if wantquick

else

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE RESAMPLING STUFF

  % define wantresampleruns and resampling
  switch options.xvalmode
  case 0
    wantresampleruns = [];
    resampling = 0;
  case 1
    wantresampleruns = 1;
    half1 = copymatrix(zeros(1,length(data)),1:round(length(data)/2),1);
    half2 = ~half1;
    resampling = [(1)*half1 + (-1)*half2;
                  (-1)*half1 + (1)*half2];
  case 2
    wantresampleruns = 0;
    resampling = [];
    for p=1:length(data)
      half1 = copymatrix(zeros(1,size(data{p},2)),1:round(size(data{p},2)/2),1);
      half2 = ~half1;
      resampling = cat(2,resampling,[(1)*half1 + (-1)*half2;
                                     (-1)*half1 + (1)*half2]);
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE STIMULUS AND DATA

  %%%%% CLUSTER CASE

  if usecluster

    % save stimulus
    while 1
      localfile0 = [scratchdir '/' sprintf('stim%d_%s.mat',subject,randomword(5))];
      if ~exist(localfile0,'file')
        break;
      end
    end
    save(localfile0,'stimulus');
    filestodelete{end+1} = localfile0;
    clear stimulus;
    
    % define stimulus
    stimulus = @() loadmulti(localfile0,'stimulus');

    % save data
    while 1
      % directory name that will contain 001.bin, etc.
      localfile0 = [scratchdir '/' sprintf('data%d_%s',subject,randomword(5))];
      if ~exist(localfile0,'dir')
        break;
      end
    end
    assert(mkdir(localfile0));
    for p=1:numruns
      savebinary([localfile0 sprintf('/%03d.bin',p)],'single',squish(data{p},dimdata)');  % notice squish
    end
    filestodelete{end+1} = localfile0;
    clear data;

    % define data
    binfiles = cellfun(@(x) [localfile0 sprintf('/%03d.bin',x)],num2cell(1:numruns),'UniformOutput',0);
    data = @(vxs) cellfun(@(x) double(loadbinary(x,'single',[0 numvxs],-vxs)),binfiles,'UniformOutput',0);

    % prepare the output directory
    while 1
      localfile0 = [scratchdir '/' sprintf('prfresults%d_%s',subject,randomword(5))];
      if ~exist(localfile0,'dir')
        break;
      end
    end
    filestodelete{end+1} = localfile0;
    filestodelete{end+1} = [localfile0 '.mat'];  % after consolidation
    outputdir = localfile0;

  %%%%% NON-CLUSTER CASE

  else

    stimulus = {stimulus};
    data = @(vxs) cellfun(@(x) subscript(squish(x,dimdata),{vxs ':'})',data,'UniformOutput',0);
    outputdir = [];

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE OPTIONS

  % last-minute prep
  if iscell(noisereg)
    noiseregINPUT = {noisereg};
  else
    noiseregINPUT = noisereg;
  end

  % construct the options struct
  opt = struct( ...
    'outputdir',outputdir, ...
    'stimulus',stimulus, ...
    'data',data, ...
    'vxs',options.vxs, ...
    'model',{model}, ...
    'seed',seedfun, ...
    'optimoptions',{{'Display' options.display 'Algorithm' 'levenberg-marquardt' 'MaxIter' options.maxiter}}, ...
    'wantresampleruns',wantresampleruns, ...
    'resampling',resampling, ...
    'metric',@calccod, ...
    'maxpolydeg',options.maxpolydeg, ...
    'wantremovepoly',1, ...
    'extraregressors',noiseregINPUT, ...
    'wantremoveextra',0, ...
    'dontsave',{{'modelfit' 'opt' 'vxsfull' 'modelpred' 'testdata'}});  % 'resnorms' 'numiters' 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIT MODEL

  %%%%% CLUSTER CASE

  if usecluster

    % define job name
    jobname = sprintf('%d_%s',subject,makedirid(opt.outputdir,1));

    % submit jobs
    jobids = CHPCsubmit(userid,jobname,'run_fitnonlinearmodel.sh',options.numperjob, ...
                        1,ceil(length(options.vxs)/options.numperjob),[], ...
                        {'data' 'stimulus' 'bad' 'd' 'xx' 'yy' 'modelfun' 'model'});

    % record additional files to delete
    filestodelete{end+1} = sprintf('~/sgeoutput/job_%s.*',jobname);  % .o and .e files
    filestodelete{end+1} = sprintf('~/mcc/job_%s.mat',jobname);

    % wait for jobs to finish
    CHPCwait(jobname,jobids,userid);

    % consolidate and load the results
    fitnonlinearmodel_consolidate(outputdir);
    a1 = load([outputdir '.mat']);

  %%%%% NON-CLUSTER CASE

  else

    a1 = fitnonlinearmodel(opt);

  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE OUTPUT

% depending on which analysis we did (quick or full optimization),
% we have to get the outputs in a common format
if wantquick
  paramsA = permute(squish(supergridseeds,dimdata),[3 2 1]);  % fits x parameters x voxels
  rA = squish(rvalues,dimdata)';                              % fits x voxels
else
  paramsA = a1.params;                                        % fits x parameters x voxels
  rA = a1.trainperformance;                                   % fits x voxels
end

% calc
numfits = size(paramsA,1);

% init
clear results;
results.ang =      NaN*zeros(numvxs,numfits);
results.ecc =      NaN*zeros(numvxs,numfits);
results.expt =     NaN*zeros(numvxs,numfits);
results.rfsize =   NaN*zeros(numvxs,numfits);
results.R2 =       NaN*zeros(numvxs,numfits);
results.gain =     NaN*zeros(numvxs,numfits);
results.resnorms = cell(numvxs,1);
results.numiters = cell(numvxs,1);

% massage model parameters for output and put in 'results' struct
results.ang(options.vxs,:) =    permute(mod(atan2((1+res(1))/2 - paramsA(:,1,:), ...
                                                  paramsA(:,2,:) - (1+res(2))/2),2*pi)/pi*180,[3 1 2]);
results.ecc(options.vxs,:) =    permute(sqrt(((1+res(1))/2 - paramsA(:,1,:)).^2 + ...
                                             (paramsA(:,2,:) - (1+res(2))/2).^2),[3 1 2]);
results.expt(options.vxs,:) =   permute(posrect(paramsA(:,5,:)),[3 1 2]);
results.rfsize(options.vxs,:) = permute(abs(paramsA(:,3,:)) ./ sqrt(posrect(paramsA(:,5,:))),[3 1 2]);
results.R2(options.vxs,:) =     permute(rA,[2 1]);
results.gain(options.vxs,:) =   permute(posrect(paramsA(:,4,:)),[3 1 2]);
if ~wantquick
  results.resnorms(options.vxs) = a1.resnorms;
  results.numiters(options.vxs) = a1.numiters;
end

% reshape
results.ang =      reshape(results.ang,      [xyzsize numfits]);
results.ecc =      reshape(results.ecc,      [xyzsize numfits]);
results.expt =     reshape(results.expt,     [xyzsize numfits]);
results.rfsize =   reshape(results.rfsize,   [xyzsize numfits]);
results.R2 =       reshape(results.R2,       [xyzsize numfits]);
results.gain =     reshape(results.gain,     [xyzsize numfits]);
results.resnorms = reshape(results.resnorms, [xyzsize 1]);
results.numiters = reshape(results.numiters, [xyzsize 1]);

% add some more stuff
results.meanvol =  meanvol;
results.noisereg = noisereg;
results.params =   paramsA;
results.options = options;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLEAN UP

% no clean up necessary in the quick case
if ~wantquick

  %%%%% CLUSTER CASE

  if usecluster

    % delete local files and directories  [should make this a function!]
    for p=1:length(filestodelete)
      if exist(filestodelete{p},'dir')  % first dir, then file
        rmdir(filestodelete{p},'s');
      elseif exist(filestodelete{p},'file')
        delete(filestodelete{p});
      end
    end

  %%%%% NON-CLUSTER CASE

  else

  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REPORT

fprintf('*** analyzePRF_chpc: ended at %s (%.1f minutes). ***\n', ...
        datestr(now),etime(clock,stime)/60);
