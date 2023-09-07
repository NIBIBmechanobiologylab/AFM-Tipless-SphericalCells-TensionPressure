format long
clc

% ======================== %
% AFM Tipless Data Analysis
% Created by cparvini
% Based upon code created
% by Alex Cartagena-Rivera
% ======================== %
% This script runs behind the 
% TiplessGUI app, and cannot be
% run on its own. It will move
% through a chosen directory
% and analyze curves in turn,
% asking for user input before 
% continuing to the next curve.
% ======================== %

% This code finds the objects inside the GUI. We use these "handles" to
% reference these objects later and update them! (e.g. we can use
% headerHandle.Text = '...' to change the header on the GUI window)
figHandle = app.Plot1;
if ~isempty(figHandle)
    axesHandle = findobj(figHandle,'Type','axes');
end
cla(axesHandle)
headerHandle = app.DataLabel;

%% Check Directory for Data
% Grab the pathname from the GUI window
pathname = app.pathname;

% Search for NanoScope text files AND Bruker .spm files AND JPK text files
Files = dir([pathname '\*NanoScope*.txt']);
FilesTemp = dir([pathname '\*.*']);
brukerEndings = "." + digitsPattern(3);
FilesTemp([~endsWith({FilesTemp.name},brukerEndings)]) = [];
FilesJPK = dir([pathname '\*qi-data*.txt']);
FilesJPK2 = dir([pathname '\*force-save*.txt']);

% Combine our lists
Files = [Files;FilesTemp;FilesJPK;FilesJPK2];

% Remove directories from the analysis, if they're in the list
Files([Files.isdir]) = [];

%% Make the log file
fid = fopen([pathname sprintf('\\TiplessAnalysisLog_%s.txt',date)],'w');
fprintf(fid,'Tipless AFM Data Analysis, %s\r\n=================\r\n\r\n', date);

FileLabel = cell(length(Files),1);
SurfaceTension = cell(size(FileLabel));
ContactPoint = cell(size(FileLabel));
RSquared = cell(size(FileLabel));

%% Begin looping through the files
for i = 1:length(Files)
    
    clearvars -except app fid i Files pathname figHandle headerHandle ...
        FileLabel SurfaceTension ContactPoint RSquared lb ub
    
    % Reset the plot for the user
    figHandle = app.Plot1;
    if ~isempty(figHandle)
        axesHandle = findobj(figHandle,'Type','axes');
    end
    cla(axesHandle)
    
    % Do a little pre-processing to figure out what filetype we're handling
    filepath = [Files(i).folder '\' Files(i).name];
    [~,~,fileExt] = fileparts(filepath);
    fileExt = strrep(fileExt,'.','');
    
    % Load File Data
    switch true
        case strcmpi(fileExt,'txt')
            if contains(Files(i).name,'NanoScope')
                % .txt files extracted from Bruker files
                rawData = importdata(filepath);
                Data = rawData.data;
                z = Data(:,1).*(1e3);           % Load Z extension [nm]
                defl = Data(:,2);               % Load Deflection [nm]
            elseif contains(Files(i).name,{'qi-data','force-save'})
                % .txt files extracted from JPK files
                allData = importdata(filepath,'\n');
                headerLineCount = 0;
                idx = [];
                for j = 1:size(allData,1)
                    if contains(allData{j},'#')
                        if contains(allData{j},'columns:')
                            colNames = allData{j};
                        end
                        if contains(allData{j},'springConstant:')
                            stiffVal = allData{j};
                        end
                        headerLineCount = headerLineCount+1;
                    else
                        break;
                    end
                end
                
                rawData = importdata(filepath, ' ', headerLineCount);
                temp = strsplit(colNames,': ');
                colNamesArray = strsplit(temp{2},' ');
                
                temp = strsplit(stiffVal,': ');
                k = str2double(strsplit(temp{2},' '));
                
                Data = rawData.data;
                z = Data(:,find(strcmpi(colNamesArray,'verticalTipPosition'),1))/(1e-9);        % Load Z extension [nm]
                z = flip(z);
                defl = Data(:,find(strcmpi(colNamesArray,'vDeflection'),1))/(1e-9)/k;     % Load Deflection [nm]
                
            end
        case endsWith(fileExt,digitsPattern(3))
            % .spm Bruker files
            [z,defl,time,dt,k] = loadBrukerAFMData(filepath,app.measurementMedium);
            z = z ./ 1e-9;
            defl = defl ./ 1e-9;
            
            if abs(diff(defl(1:2))) > abs(10*median(gradient(defl(1:100))))
                z(1) = [];
                defl(1) = [];
            end
            
            if mean(z(1:100)) <= 0
                z = z - z(1) + (z(2)-z(1));
            end
            
        otherwise
            error('A file that was designated for analysis has an unknown extension!');
    end
    
    r = length(z);                  % Length of the Input Data
    if ~exist('lb')
        lb = round(0.7*r);              % Lower Bound for Fitting
        ub = round(0.9*r);              % Upper Bound for Fitting
    else
        if ub > r
            ub = r;
        end
        if lb > ub
            lb = 0.7*ub;
        end
    end
    T0 = 1e-6;                      % Initial Guess for the Tension
    
    app.EndIndexSpinner.Limits(2) = r;
    
    % Update the command window
    fprintf('Analyzing File: %s\n',Files(i).name);
    
    % Update the header for the user
    headerHandle.Text = sprintf('Current File: %d of %d',i,length(Files));
    
    % Parameters
    if isempty(app.cellRadius)
        R = input(sprintf('Please enter the Cell Radius [nm]: ')); % Spherical Cell Radius (nm)
    else
        R = app.cellRadius; % Spherical Cell Radius (nm)
    end
    if ~exist('k','var')
        k = app.stiffness; % Calibrated Cantilever Spring Constant (nN/nm)
    end
    
    % Do an initial fit
    xnew = z((lb):(ub));
    ynew = defl((lb):(ub));

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
    
    % Show us what our initial guess looks like
    figHandle = app.Plot1;
    if ~isempty(figHandle)
        axesHandle = findobj(figHandle,'Type','axes');
    end
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
    
    % Enable the inputs in the GUI for the user
    app.StartIndexSpinner.Enable = true;
    app.EndIndexSpinner.Enable = true;
    app.AcceptButton.Enable = true;
    app.RePlotButton.Enable = true;
    app.EndIndexSpinner.Value = ub;
    app.StartIndexSpinner.Value = lb;
    
    % Send data to the GUI so it can be used for re-plotting
    app.xdata = z;
    app.ydata = defl;
    app.R = R;
    app.k = k;
    
    % Hang out here while the user re-plots and fine tunes the data used
    % for fitting.   
    waitfor(app,'acceptButtonPress',true);
    
    % Reset the acceptButtonPress value before moving on
    app.acceptButtonPress = false;
    
    % Grab the user-optimized lb/ub
    lb = app.StartIndexSpinner.Value;           % Lower Bound for Fitting
    ub = app.EndIndexSpinner.Value;           % Upper Bound for Fitting
    
    % Disable the inputs in the GUI for the user
    app.StartIndexSpinner.Enable = false;
    app.EndIndexSpinner.Enable = false;
    app.AcceptButton.Enable = false;
    app.RePlotButton.Enable = false;
    
    % Do the analysis again, using the new bounds
    xnew = z((lb):(ub));
    ynew = defl((lb):(ub));
    
    % Add a header for this file
    fprintf(fid,'File %d of %d: %s\r\n', i, length(Files), Files(i).name);
    
    xnew1 = xnew - ynew; % New x because equation has y and fit algorithm does not like that
    s = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0],'Upper',[1,1e4],'Startpoint',[1 1]); %Options to use nonlinear model
    f = fittype('(a.*(x-b)).*(1+(5/(4*6800)).*(x-b))','coefficients', {'a','b'},'options',s); %Equation used to fit data
    [fC, gof] = fit(xnew1,ynew,f,'Robust','on','Display','off','TolFun',1e-12,'TolX',1e-12,'MaxIter',400,'MaxFunEvals',600); %Fit function
    
    Tnew = [fC.a]*(k/pi); % Surface tension (N/m)
    fprintf(fid,'The cell surface tension T [pN/um] = %1.3g\r\n',Tnew*(1e6));

    Znew = [fC.b]; % Contact point
    fprintf(fid,'The contact point Z0 = %1.3g\r\n',Znew);
    
    rsq = gof.adjrsquare; %Compute adjusted R-squared value for the fit
    fprintf(fid,'The adjusted R-squared R2 = %1.3f\r\n\r\n',rsq);
    
    % Save data into table columns
    temp = strsplit(Files(i).name,{'_','.'});
    idx1 = find(contains(temp,'Cell'));
    idx2 = find(contains(temp,'Dish'));
    
    switch true
        case strcmpi(fileExt,'txt')
            if contains(Files(i).name,'NanoScope')
                % .txt files extracted from Bruker files
                temp2 = strfind(Files(i).name,digitsPattern(3));
                idx3 = temp2(1);
                FileLabel{i,1} = [temp{idx2} '-' temp{idx1} '-' Files(i).name(idx3:(idx3+2))];
            elseif contains(Files(i).name,{'qi-data','force-save'})
                % .txt files extracted from JPK files
                FileLabel{i,1} = [temp{idx2} '-' temp{idx1} '-' num2str(i)];
            end
        case endsWith(fileExt,digitsPattern(3))
            % .spm Bruker files
            idx3 = length(temp);
            FileLabel{i,1} = [temp{idx2} '-' temp{idx1} '-' temp{idx3}];
    end
    
    SurfaceTension{i,1} = Tnew*(1e6);
    ContactPoint{i,1} = Znew;
    RSquared{i,1} = rsq;
    
end

% Update the header for the user
headerHandle.Text = 'Analysis Complete!';
fclose(fid);

% Make the table and save it
resultsTable = table(FileLabel,SurfaceTension,ContactPoint,RSquared);
outputFilename = [pathname sprintf('\\TiplessAnalysisLog_%s.xlsx',date)];
writetable(resultsTable,outputFilename);

% Reset the plot for the user
figHandle = app.Plot1;
if ~isempty(figHandle)
    axesHandle = findobj(figHandle,'Type','axes');
end
cla(axesHandle)

% Open the analysis directory
winopen(pathname);