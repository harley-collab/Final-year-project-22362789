function ring_sim_4_PCB()
 tic;
  basePlotDir = 'plots_4_3_1';
  
  dataDir     = 'data_4_3_1';

  if ~exist(basePlotDir,'dir'), mkdir(basePlotDir); end
  if ~exist(dataDir,'dir'), mkdir(dataDir); end


p = gcp('nocreate');
if isempty(p)
    parpool('local', feature('numcores'));
end

C1_seg1 = linspace(0.001, 0.01, 20);%drag pairwise

C1_seg2 = linspace(0.01, 0.1, 20);

C1_seg3 = linspace(0.1, 1, 20);

C1_values = [C1_seg1, C1_seg2(2:end), C1_seg3(2:end)];%merge no overlap 

C2_seg1 = linspace(0.001, 0.01, 100);

C2_seg2 = linspace(0.01, 0.1, 100);

C2_seg3 = linspace(0.1, 1, 100);

C2_values = [C2_seg1, C2_seg2(2:end), C2_seg3(2:end)];

  m  = 1.0;
  R  = 1.025;
  h  = 0.88;
  g  = 9.81;
  I1 = 0.5772;
  I3 = 1.0253;

  r_Gc = [0; -h/2; R];

  theta0   = 0.55;
  phi_dot0 = 10;
  psi_dot0 = 0.2;

  tspan    = [0 600];

  w10 = 0;
  w20 = psi_dot0 + phi_dot0*sin(theta0);
  w30 = phi_dot0*cos(theta0);
  y0 = [0; 0; R*cos(theta0); w10; w20; w30; theta0; 0; 0];

q = parallel.pool.DataQueue;
completed = 0;
afterEach(q, @(~) progressUpdate());

function progressUpdate()
    completed = completed + 1;
    fprintf('Progress: %d / %d  (%.1f%%)\n', completed, nJobs, 100*completed/nJobs);
end
  
[iC1_grid, iC2_grid] = ndgrid(1:length(C1_values), 1:length(C2_values));
nJobs = numel(iC1_grid);
C3_crit_flat = nan(nJobs, 1);

parfor job = 1:nJobs
    iC1 = iC1_grid(job);
    iC2 = iC2_grid(job);
    C1 = C1_values(iC1);
    C2 = C2_values(iC2);

    C3_crit = findCriticalC3_NR(theta0, C1, C2, m, R, h, g, I1, I3, phi_dot0, psi_dot0, tspan, r_Gc, y0);

    C3_crit_flat(job) = C3_crit;
    send(q, true);

    C3 = C3_crit;

    options = odeset('RelTol',1e-5,'AbsTol',1e-5,'MaxStep',1e-1, 'Events', @(t,y) combinedEvents(t,y));
    [t, y] = ode45(@(t,y) eom(t,y,m,R,h,g,I1,I3,C1,C2,C3,r_Gc), tspan, y0, options);

    Ff = zeros(size(t));
    for k = 1:length(t)
        Ff(k) = turntrack(y(k,:),m,R,h,g,I1,I3,C1,C2,C3,r_Gc);
    end

    tag = sprintf('C1 %.3f C2 %.3f C3 %.10f', C1, C2, C3);

    
    dataFile = fullfile(dataDir, [tag '.dat']);
    fid = fopen(dataFile,'w');
    fprintf(fid,'# Critical C3 sweep\n');
    fprintf(fid,'# C1=%.6f  C2=%.6f  Critical C3=%.10f\n', C1, C2, C3_crit);
    fprintf(fid,'# t  rcx rcy rcz  w1 w2 w3  theta phi psi  Ff\n');
    fclose(fid);
    writematrix([t y Ff], dataFile, 'Delimiter',' ', 'WriteMode','append');
end


C3_crit_matrix = reshape(C3_crit_flat, length(C1_values), length(C2_values));

matrixFile = fullfile(dataDir, 'C3_crit_matrix_top.dat');%top PCB
%matrixFile = fullfile(dataDir, 'C3_crit_matrix_bottom.dat');%bottom PCB
fid = fopen(matrixFile, 'w');
fprintf(fid, '# C3 critical matrix\n');
fprintf(fid, '# C1  C2  C3_crit\n');
for iC1 = 1:length(C1_values)
    for iC2 = 1:length(C2_values)
        fprintf(fid, '%.6f %.6f %.10f\n', C1_values(iC1), C2_values(iC2), C3_crit_matrix(iC1,iC2));
    end
end
fclose(fid);
  toc;
end