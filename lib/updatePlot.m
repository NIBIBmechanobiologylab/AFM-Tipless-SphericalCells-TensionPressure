function [] = updatePlot(figHandle,z,defl,lb,ub,R,k)
%UPDATEPLOT Update the GUI plot for tiplessGUI
%   Update the GUI plot for tiplessGUI
    
    if ~isempty(figHandle)
        axesHandle = findobj(figHandle,'Type','axes');
    end
    
    % Grab our data in the new bounds
    xnew = z((lb):(ub));
    ynew = defl((lb):(ub));
    
    % Analyze the data between those bounds
    d0 = mean(ynew(1:round(numel(ynew)*0.4))); %Initial guess for y-axis
    xnew1 = xnew - ynew; % New x because equation has y and fit algorithm does not like that
    s = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0],'Upper',[1,1e4],'Startpoint',[1 1]); %Options to use nonlinear model
    f = fittype('(a.*(x-b)).*(1+(5/(4*6800)).*(x-b))','coefficients', {'a','b'},'options',s); %Equation used to fit data
    [fC, gof] = fit(xnew1,ynew,f,'Robust','on','Display','off','TolFun',1e-12,'TolX',1e-12,'MaxIter',400,'MaxFunEvals',600); %Fit function
    Tnew = [fC.a]*(k/pi); % Surface tension (N/m)
    Znew = [fC.b]; % Contact point
    rsq = gof.adjrsquare; %Compute adjusted R-squared value for the fit
    yfit = ((pi*Tnew)/k).*(xnew1-Znew).*(1+(5/4).*(xnew1-Znew)./(R));
    T = Tnew;
    z2 = (z-Znew)*10^(-3);
    yexact = exact(z2,R,k,T);
    
    % Update the plot
    cla(axesHandle)
    hold(axesHandle, 'on')
    plot(axesHandle, z, defl, 'b')
    dtip = plot(axesHandle, xnew, ynew, 'o');
    plot(axesHandle, xnew, yfit, 'r', 'linewidth', 5)
    plot(axesHandle, z, yexact, 'g', 'linewidth', 3)
    xlabel(axesHandle, 'Piezo extension Z (nm)')
    ylabel(axesHandle, 'Cantilever deflection d (nm)')
    title(axesHandle, 'Deflection vs. Z Extension')
    legend(axesHandle, 'Data', 'Selection', 'Fit', 'Exact', 'location', 'southeast')
    text(axesHandle, 0.05, 0.85, sprintf('$$R^2 = %1.3f$$\n$$\\Delta Z = %3.1f nm$$',...
        rsq, xnew(end)-xnew(1)),...
        'units', 'normalized','interpreter','latex','fontsize',24)
    hold(axesHandle, 'off')
    
end