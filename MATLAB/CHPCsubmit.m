function jobids = CHPCsubmit(userid,jobname,scriptname,numvxs,beginix,endix,flags,varstoexclude)

% function jobids = CHPCsubmit(userid,jobname,scriptname,numvxs,beginix,endix,flags,varstoexclude)
% 
% <userid> is a string. this should be the user who is submitting the jobs.
% <jobname> is a string with the name of the job
% <scriptname> is the name of the script (e.g. 'run_fitnonlinearmodel.sh')
% <numvxs> is the number of voxels to process in each job
% <beginix> is starting voxel index
% <endix> is ending voxel index
% <flags> (optional) is a string with some qsub flags.
%   Default: '-l nodes=1:ppn=1,walltime=4:00:00,vmem=2gb -q pe1_iD_4hr'.
% <varstoexclude> (optional) is a cell vector of variables to exclude 
%   when saving to the temporary .mat file
%
% submit jobs and return a vector of job ids.
% 
% some special things we do:
% (1) we make available the caller's workspace to the jobs
% (2) we submit in chunks, being careful not to overload the queue

% internal constants
maxsize = 10;    % warn if workspace is bigger than this number of MB
maxjobs = 1500;  % max jobs to have running at one time
jobchunk = 1000; % how many jobs to submit in each iteration
pausetime = 10;  % seconds to pause when unexpected exit code
overloadtime = 60;  % seconds to wait when overloaded before trying again

% inputs
if ~exist('flags','var') || isempty(flags)
  flags = '-l nodes=1:ppn=1,walltime=4:00:00,vmem=2gb -q pe1_iD_4hr';  % -q dque
end
if ~exist('varstoexclude','var') || isempty(varstoexclude)
  varstoexclude = {};
end

% setup
mkdirquiet('~/sgeoutput');
mkdirquiet('~/mcc');
%%%%% ensure uniqueness

if exist(sprintf('~/mcc/job_%s.mat',jobname),'file')
  error(sprintf('jobname %s already exists!',jobname));
end

%%%%% save the caller's workspace to the special file

% check size of workspace
a = evalin('caller','whos');
ok = ~ismember(cat(2,{a.name}),varstoexclude);
a = a(ok);
workspacesize = sum(cat(1,a.bytes));

% give warning if the workspace is big
if workspacesize/1000/1000 > maxsize
  warning(sprintf('We are saving to disk a workspace that is larger than %.1f MB!',maxsize));
end

% save caller's workspace to the special .mat file
fprintf('saving workspace to disk...');
evalin('caller',sprintf('saveexcept(''~/mcc/job_%s.mat'',%s);',jobname,cell2str(varstoexclude)));
fprintf('done.\n');

% do in chunks
jobids = [];
totalchunks = ceil((endix-beginix+1)/jobchunk);
for zz=1:totalchunks
  fprintf('working on chunk %d of %d.\n',zz,totalchunks);

  %%%%% wait until not overloaded

  fprintf('checking to see if CHPC is overloaded.\n');
  while 1
    [s,r] = unix(sprintf('qstat -t -u %s | wc -l',userid));
    fprintf('exit status was %d.\n',s);
    while s~=0
      fprintf('STATUS WASNT 0, WEIRD; TRY AGAIN IN A FEW SECONDS\n');
      pause(pausetime);
      [s,r] = unix(sprintf('qstat -t -u %s | wc -l',userid));
      fprintf('exit status was %d.\n',s);
    end
    if str2double(r) > maxjobs
      fprintf('overloaded, so waiting (%s).\n',datestr(now));
      pause(overloadtime);
    else
      break;
    end
  end
  fprintf('CHPC is not overloaded right now, so let''s proceed.\n');

  %%%%% issue the qsub call on chpc

  beginix0 = (zz-1)*jobchunk + beginix;
  endix0 = min(beginix0 + jobchunk - 1,endix);
  cmd = sprintf('qsub -t %d-%d -v MYSCRIPT="%s",MYARG1="~/mcc/job_%s.mat",MYARG2="%d" -N job_%s -o ~/sgeoutput/job_%s.o -e ~/sgeoutput/job_%s.e %s ~/mcc/matlabsge.sh', ...
    beginix0,endix0,scriptname,jobname,numvxs,jobname,jobname,jobname,flags);
  fprintf('this is the qsub command:\n\n%s\n\n',cmd);
  fprintf('issuing the qsub command on chpc...');
  [s,r] = unix(cmd);
  fprintf('exit status was %d.\n',s);
  while s~=0
    fprintf('STATUS WASNT 0, WEIRD; TRY AGAIN IN A FEW SECONDS\n');
    pause(pausetime);
    [s,r] = unix(cmd);
  end
  jobids = [jobids str2double(regexp(r,'\d+','match'))];
  fprintf('done.\n');

  %%%%% wait to take effect

  if zz ~= totalchunks
    fprintf('waiting before trying next chunk.\n');
    pause(overloadtime);
  end

end
