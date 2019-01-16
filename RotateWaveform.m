clear all;

% read in bundle
dir = '/nethome/numrel/datafiles/Waveforms/LL-series/q6_LL_D9_a0.6_th1_45_th2_225_m320';
bn = modeBundle( dir, 2:3, 75 );

% read in kick
file = sprintf('%s/kicksr75.asc', dir );
data = dlmread(file,'\t');
time = data(end,1);
kick = data(end,2:4);

% setup mesh grid
[TH, PH] = skyMath.setup_mesh_grid( pi/20, 2*pi );

[rs, cs] = size(TH);

count = 1;
tcount = cs;
for j=1
   for k=1:cs
      fprintf('%d/%d...\n',count,tcount); 
      % setup frame
      frm = frame( 0, PH(j,k), TH(j,k), 0, 'tmp' );

      % rotate by frame: ie take the z-axis and place it at angles 
      % defined by frm
      ro = bn.rotate( frm );
      
      %sm  = bn.recomposeMeshAtSingleTime( 200, pi/60, pi/60);
      %sm2 = ro.recomposeMeshAtSingleTime( 200, pi/60, pi/60);

      % rotate kick vector.  since we are using a fixed frame, we 
      % want the inverse rotation to be consistent
      k_rot = frm.rotate_vector_to_frame( time, kick );

      ampl(count) = max(ro.Ampl(:,5))/max(ro.Ampl(:,1));
      kz  (count) = k_rot(3);
      count       = count + 1;
   end
end

figure(1); clf;
semilogy( kz, ampl, '*b' );

    
