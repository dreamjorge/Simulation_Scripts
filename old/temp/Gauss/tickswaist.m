function tickswaist(nrticksy,window,so,wo)

  stepticks  = window/(nrticksy-1);
  stepso     = so/window;
  f =@waistGaussianBeam;
  %nrticksy     = 2*timesso+1;                                                   % Number of ticks symetric of 0.   
  yticksv      = zeros(1,nrticksy);                                           % Vector for values in ticks.
  yticklabelsv = cell(1,nrticksy);
ii=1;
  for jj = floor(nrticksy/2)+3:1:nrticksy
      mm                                   = jj - floor(nrticksy/2)-1;
      ntimesso                             = mm;
      yticksv(jj)                          = f(ii*stepticks*so,wo.^2);                                             % How many times of so evalued in Waist.
      yticklabelsv{jj}                     = ['$$',latex(simplify(f(ii*stepticks*sym('wo').^2,sym('wo')))),'$$'];
      yticksv(floor(nrticksy/2)+1-mm)      = -yticksv(jj); 
      yticklabelsv{floor(nrticksy/2)+1-mm} = ['$$-',latex(simplify(f(ii*stepticks*sym('wo').^2,sym('wo')))),'$$'];
      ii=ii+1;
  end

  yticksv(floor(nrticksy/2)+2)  = f(0,wo.^2);
  yticklabelsv{floor(nrticksy/2)+2}   = '$wo$';
  yticksv(floor(nrticksy/2))    = -f(0,wo.^2);
   yticklabelsv{floor(nrticksy/2)}   = '$-wo$';
   yticklabelsv{floor(nrticksy/2)+1} = '$$0$$';
  %%
  set(gca,'ytick',yticksv);                                                   % Set values of ticks. 
  set(gca,'yticklabel',yticklabelsv)  
end