function ticksx(nrticksx,window,so)
  stepticks = window/(nrticksx-1);                                      % Number of ticks symetric of 0.
  xticksv      = zeros(1,nrticksx);                                           % Vector for values in ticks.
  xticklabelsv = cell(1,nrticksx);                                            % Cell for strings of labels in ticks.
  %cycle for how many times of so
  for ii = 1:nrticksx
      xticksv(ii)      = (ii-floor(nrticksx/2)-1)*stepticks;                         % How many times of so.
      if ( xticksv(ii) == 0)                                                  % If x is zero only put 0 in xlabel.
           xticklabelsv{ii} = 0;                                              
      else                                                                    % Else It takes value of times so.
           xticklabelsv{ii} = [num2str(xticksv(ii)),'$s_0$'];                    
      end
  end
  set(gca,'xtick',xticksv);                                                   % Set values of ticks. 
  set(gca,'xticklabel',xticklabelsv)   