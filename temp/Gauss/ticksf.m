function ticksf(f,xstep,xunit,nrticksx,nrticksy) 

%generate x-ticks in terms of units of x                                           % Number of ticks symetric of 0.
xticksv      = zeros(1,nrticksx);                                           % Vector for values in ticks.
xticklabelsv = cell(1,nrticksx);                                            % Cell for strings of labels in ticks.

%cycle for how many times of so
for ii = 1:nrticksx
    xticksv(ii)      = (ii-floor(nrticksx/2)-1)*xstep;                      % How many times of so.
    if ( xticksv(ii) == 0)                                                  % If x is zero only put 0 in xlabel.
         xticklabelsv{ii} = 0;                                              
    else                                                                    % Else It takes value of times so.
         xticklabelsv{ii} = [num2str(xticksv(ii)),xunit];                    
    end
end

set(gca,'xtick',xticksv);                                                   % Set values of ticks. 
set(gca,'xticklabel',xticklabelsv)                                          % Set labels in ticks. 

%generate y-ticks in terms of y-units                                       % Number of ticks symetric of 0.   
yticksv      = zeros(1,nrticksx);                                           % Vector for values in ticks.
yticklabelsv = cell(1,nrticksx);     
% Cell for strings of labels in ticks.

    
%cycle for how many times of so is evaluated in Radius functions
for jj =1:nrticksy
     % How many times of so evalued in Radius.
    
  if (strcmp(xunit,'so'))

    fvalue      = f((jj-floor(nrticksx/2)-1)*xstep,1);
    fticks      = f((jj-floor(nrticksx/2)-1)*xstep,1);
    yticksv(jj) = fvalue;
    if isnumeric(fticks)
        if fticks == inf
           fticks = 0; 
        end
      fticks = [rats(fticks),'so'];
    else
      fticks = [latex(fticks),'so'];
    end 

  elseif(strcmp(xunit,'wo'))
  
  
  end
     
    
    yticklabelsv{jj} = ['$$',fticks,'$$']          

end

yticksv(yticksv==inf)=0;
set(gca,'ytick',yticksv);                                                  
set(gca,'yticklabel',yticklabelsv)                                          % Set labels in ticks. 

end