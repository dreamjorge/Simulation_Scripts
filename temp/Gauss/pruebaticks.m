function tickswaist(nrticksy)

  f =@waistGaussianBeam;
  %nrticksy     = 2*timesso+1;                                                   % Number of ticks symetric of 0.   
  yticksv      = zeros(1,nrticksy);                                           % Vector for values in ticks.
  yticklabelsv = cell(1,nrticksy);

  for jj = floor(nrticksy/2)+2:1:nrticksy
      mm                                   = jj - floor(nrticksy/2)-1;
      ntimesso                             = mm;
      yticksv(jj)                          = f((ntimesso-1)*so,wo.^2);                                             % How many times of so evalued in Waist.
      yticklabelsv{jj}                     = ['$$',latex(simplify(f((ntimesso-1)*sym('wo').^2,sym('wo')))),'$$'];
      yticksv(floor(nrticksy/2)+1-mm)      = -yticksv(jj); 
      yticklabelsv{floor(nrticksy/2)+1-mm} = ['$$-',latex(simplify(f(-(ntimesso-1)*sym('wo').^2,sym('wo')))),'$$'];

  end

   yticklabelsv{floor(nrticksy/2)+1} = '$$0$$';
  %%
  set(gca,'ytick',yticksv);                                                   % Set values of ticks. 
  set(gca,'yticklabel',yticklabelsv)   