function run_all(what)
%RUN_ALL  Single entry point for every experiment in this repository.
%
%   run_all              run the self-check and all figures (about 12 min)
%   run_all('selftest')  kernel self-check, 4 invariants            ~5 s
%   run_all('fig2a')     Fig. 2(a): gain vs bandwidth split, 2 UEs   ~40 s
%   run_all('fig2b')     Fig. 2(b): 5 UEs, per-user gain            ~40 s
%   run_all('fig3')      Fig. 3 + Fig. 4: patterns and heat maps    ~50 s
%   run_all('fig5')      Fig. 5: greedy (Alg. 1) and Adam GD         ~3.6 min
%   run_all('fig5conv')  only the equivalence and start-point checks ~1.5 min
%   run_all('fig6')      Fig. 6 + Fig. 7: comparison against [5]     ~2.0 min
%
%   Every job writes its figures and tables into matlab/out/.

if nargin < 1 || isempty(what), what = 'all'; end
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here,'core'));
addpath(fullfile(here,'alg'));

% key, script, kind ('f' = function file, 's' = script file), description
jobs = {
    'selftest', 'run_p0',     'f', 'kernel self-check: 4 invariants (P0)'
    'fig2a',    'run_a1',     'f', 'Fig. 2(a): gain vs bandwidth split, 2 UEs'
    'fig2b',    'run_a3',     'f', 'Fig. 2(b): 5 UEs, per-user gain and fairness'
    'fig3',     'run_a45',    'f', 'Fig. 3 (azimuth-elevation pattern) and Fig. 4 (azimuth-frequency)'
    'fig5',     'run_a6',     'f', 'Fig. 5: greedy Algorithm 1 and Adam gradient descent'
    'fig5conv', 'run_a6_conv','f', 'Fig. 5 support: equivalence check and start-point diagnostic'
    'fig6',     'run_a7',     'f', 'Fig. 6 and Fig. 7: comparison against Algorithm 1 of [5]'
};

if strcmpi(what,'all')
    sel = true(size(jobs,1),1);
else
    sel = strcmpi(jobs(:,1), what);
    if ~any(sel)
        fprintf('known jobs: all, %s\n', strjoin(jobs(:,1).', ', '));
        return
    end
end

for k = find(sel).'
    fprintf('\n########## %s : %s ##########\n', jobs{k,1}, jobs{k,4});
    t0 = tic;
    switch jobs{k,3}
        case 'f', feval(jobs{k,2});
        case 's', evalin('base', sprintf('run(''%s'')', fullfile(here,[jobs{k,2} '.m'])));
        %  A script file runs in the workspace of its caller, so it would overwrite the
        %  loop index and the timer of this function.  Evaluating it in the base workspace
        %  keeps run_all itself intact; all paths are global, so nothing else is needed.
    end
    fprintf('########## %s done in %.1f s ##########\n', jobs{k,1}, toc(t0));
end
