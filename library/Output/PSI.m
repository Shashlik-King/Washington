function PSI(pile,data,element,scour,settings,t,z,p,y,Q,zQ,plug_unplug,iii)
%% FILE FOR GENERATING PSI-INPUT FOR SACS
% BASED ON MANUAL FOR SACS
%--------------------------------------------------------------------------
% CHANGE LOG
% 2013.07.24        MUOE        Reprogramming to include generic project 
%                               info and t-z curves. Incorporating into 
%                               CMAT main.m for easier computation of soil 
%                               response and generation of SACS-input
%--------------------------------------------------------------------------
%% Input parameters
%--------------------------------------------------------------------------
% data.xx:          data regarding project, location, creator of data
% t.top:            t-data (for t-z curves) for top node of each element
% t.bottom:         t-data (for t-z curves) for bottom node of each element
% z:                z-data (for t-z curves) for each element (same for all
%                   elements, same for top and bottom node)
% p.top:            p-data (for p-y curves) for top node of each element
% p.bottom:         p-data (for p-y curves) for bottom node of each element
% y:                y-data (for p-y curves)
% element.xx        data regarding depth of top and bottom node for each
%                   element
%--------------------------------------------------------------------------
%% Generate output for SACS
%--------------------------------------------------------------------------
% Choose whether t-z and p-y curves should vary linearly between two depths
% or every curve should be defined as constant within a specified interval
lin_con = 2; % 1 = linear, 2 = constant
%--------------------------------------------------------------------------
%% Basic information
%--------------------------------------------------------------------------
name = ['psiinp.',data.location];
F = fopen(name,'Wt');
fprintf(F,'%s %s %s %s','*SOIL PROFILE FOR LOCATION:', data.location ,'at', data.project ,'');
fprintf(F,'\n');
fprintf(F,'********************************************************************************\n');
fprintf(F,'%s %0.2f', '*LOCAL SCOUR [m]:',scour.local);
fprintf(F,'\n');
fprintf(F,'%s %0.2f', '*WATER DEPTH [m]:',scour.water_depth);
fprintf(F,'\n');
fprintf(F,'%s %.2f', '*PILE LENGTH [m]:',pile.length(iii));
fprintf(F,'\n');
fprintf(F,'%s %.2f', '*PILE DIAMETER [m]:',pile.diameter(iii));
fprintf(F,'\n');
fprintf(F,'%s %.2f', '*PILE HEAD [m VREF]:',pile.head);
fprintf(F,'\n');
fprintf(F,'%s %.2f', '*PILE STICK-UP [m]:',pile.stick_up);
fprintf(F,'\n');
fprintf(F,'%s %s', '*FILE CREATED BY:',data.prepared_by,',',datestr(now));
%--------------------------------------------------------------------------
%% Soil layout
%--------------------------------------------------------------------------
fprintf(F,'\n********************************************************************************\n');
fprintf(F,'SOIL');

if settings.axial_loading
%--------------------------------------------------------------------------
%% t-z curves
%--------------------------------------------------------------------------
% FORMAT
% Putting t.top, t.bottom and z into two matrices, tz.top and tz.bottom,
% where the format is [t1.1 z1.1  t1.2 z1.2  t1.3 z1.3  t1.4 z1.4  t1.5 z1.5 .. ;
%                      t2.1 z2.1  t2.2 z2.2  t2.3 z2.3  t2.4 z2.4  t2.5 z2.5 .. ; ..]
% ".." resembles the possibility that there can be more than more than five
% points per depth interval and more than two depth intervals
i = 1;
m = 1;
while i<=size(t.top,1)
    for j=1:size(z.top,2);
        % Top nodes
        tz.top(m,j*2-1)     = t.top(i,j)/1E4; % to go from [kPa] to [kN/cm^2]
        tz.top(m,j*2)       = z.top(i,j)*1E2; % to go from [m] to [cm]
        % Bottom nodes
        tz.bottom(m,j*2-1)  = t.bottom(i,j)/1E4; % to go from [kPa] to [kN/cm^2]
        tz.bottom(m,j*2)    = z.bottom(i,j)*1E2; % to go from [m] to [cm]
    end
    m = m + 1;
    i = i + 1;
end
% PRINT
fprintf(F,'\n********************************************************************************\n');
fprintf(F,'*T-Z CURVES');
fprintf(F,'\n********************************************************************************\n');
fprintf(F,'%s%3.0f                    %s','SOIL TZAXIAL HEAD',2*length(element.level),'TYP1'); % t-z curves are defined in top and bottom node of each element
fprintf(F,'\n');
for i = 1:size(element.level,1)-1 % loop to print t-z curves
    for top_bottom = 1:2
        if top_bottom == 1
            tzdata = tz.top;
        else
            tzdata = tz.bottom;
        end
        n_max = 5; % maximum number of t-z points on one input line
        fprintf(F,'%s %s %s','*LAYER CODE:',element.model_axial{i},'for axial resistance and for t-z'); % Text to identify soil layer
        fprintf(F,'\n');
        fprintf(F,'%s %2.f','SOIL         SLOC   ',size(z.top,2)); % Text to inform SACS that the following line contains user defined t-z data. 'SLOCSM' if t-z curves are symmetric (SM = symmetric). 
        if size(z.top,2) > 30 || size(z.bottom,2) > 30
            error('For use of more than 30 points on the t-z curve, you must specify the number in the header')
        end
        if top_bottom == 1
            % the print of the depth interval for each set of t-z values
            % can cause problems if depth is greater than 100 m - this is
            % because SACS only reads from certain columns, and three
            % digits before the comma, e.g. 100.50, will not fit inside the
            % before mentioned columns
            if lin_con == 2 % curves are constant in the specified interval
                fprintf(F,'  %5.2f %5.2f', abs(element.level(i,1)-pile.head),abs(element.level(i,1)+(element.level(i,2)-element.level(i,1))/2-pile.head)); % Upper depth limit for t-z curve
            else % if not constant, then always varying linearly
                fprintf(F,'  %5.2f', abs(element.level(i,1)-pile.head)); % Upper depth limit for t-z curve
            end
        else
            if lin_con == 2 % curves are constant in the specified interval
                fprintf(F,'  %5.2f %5.2f', abs(element.level(i,1)+(element.level(i,2)-element.level(i,1))/2-pile.head),abs(element.level(i,2)-pile.head)); % Lower depth limit for t-z curve
            else % if not constant, then always varying linearly
                fprintf(F,'  %5.2f', abs(element.level(i,2)-pile.head)); % lower depth limit for t-z curve
            end
        end
        fprintf(F,'\n');
        for ii = 1:size(tzdata,1)
            for jj = 1:size(tzdata,2)
                if floor(jj/2) == jj/2 % even number, is valid for all z-values (lower precision is needed here compared to t-values)
                    if tzdata(ii,jj) > 0 % valid for positive numbers and 0 (need to distinguish because of "-"-sign       
                        if tzdata(ii,jj) >= 10 % number is bigger than 10 (meaning there is room for one less digit after the comma)
                            indicator(ii,jj) = 1; 
                        else % number is smaller than 10
                            indicator(ii,jj) = 2;
                        end
                    else % valid for negative numbers
                        if tzdata(ii,jj) <= -10 % number is smaller than -10
                            indicator(ii,jj) = 3;
                        else % number is bigger than -10
                            indicator(ii,jj) = 4;
                        end
                    end
                else % uneven number, is valid for all t-values (see comments further up)
                    if tzdata(ii,jj) > 0 % valid for positive numbers and 0 (need to distinguish because of "-"-sign       
                        if tzdata(ii,jj) >= 10
                            indicator(ii,jj) = 5; 
                        else
                            indicator(ii,jj) = 6;
                        end
                    else % valid for negative numbers
                        if tzdata(ii,jj) <= -10
                            indicator(ii,jj) = 7;
                        else
                            indicator(ii,jj) = 8;
                        end
                    end
                end
            end
        end
        form = {'%4.1f ';'%4.2f ';'%4.0f ';'%4.1f ';'%6.3f ';'%6.4f ';'%6.2f ';'%6.3f '}; % form indicator is connected to numbers 1-4 shown above
        for nsec = 1:ceil(size(tzdata,2)/(2*n_max)); % number of sections to print, maximum n_max t-z pairs in each section
            fprintf(F,'SOIL         T-Z ');
            form_indicator = [form{indicator(i,1+2*n_max*(nsec-1):min(2*n_max+2*n_max*(nsec-1),size(tzdata,2)))}]; % chooses format for each line of t-z value pairs
            fprintf(F,form_indicator, tzdata(i,1+2*n_max*(nsec-1):min(2*n_max+2*n_max*(nsec-1),size(tzdata,2)))); % prints the t-z value pairs
            fprintf(F,'\n');
        end
    end
end
%--------------------------------------------------------------------------
%% Q-z curves
%--------------------------------------------------------------------------
% FORMAT
% Putting Q.top, Q.bottom and zQ into two matrices, tz.top and tz.bottom,
% where the format is [Q1.1 zQ1.1  Q1.2 zQ1.2  Q1.3 zQ1.3  Q1.4 zQ1.4  Q1.5 zQ1.5 .. ;
%                      Q2.1 zQ2.1  Q2.2 zQ2.2  Q2.3 zQ2.3  Q2.4 zQ2,4  Q2.5 zQ2.5 .. ; ..]
% ".." resembles the possibility that there can be more than more than five
% points per depth interval and more than two depth intervals
i = 1;
m = 1;
while i<=size(Q.top,1)
    for j=1:size(zQ.top,2);
        % Top nodes
        Qz.top(m,j*2-1)     = Q.top(i,j)/1E4; % to go from [kPa] to [kN/cm^2]
        Qz.top(m,j*2)       = zQ.top(i,j)*1E2; % to go from [m] to [cm]
        % Bottom nodes
        Qz.bottom(m,j*2-1)  = Q.bottom(i,j)/1E4; % to go from [kPa] to [kN/cm^2]
        Qz.bottom(m,j*2)    = zQ.bottom(i,j)*1E2; % to go from [m] to [cm]
    end
    m = m + 1;
    i = i + 1;
end
% PRINT
fprintf(F,'********************************************************************************\n');
fprintf(F,'*Q-Z CURVES');
fprintf(F,'\n********************************************************************************\n');
fprintf(F,'%s%3.0f                    %s','SOIL BEARING HEAD',1,'TYP1');
fprintf(F,'\n');
for i = size(element.level,1)-1 % loop to print Q-z curves
    Qzdata = Qz.bottom;
    n_max = 5; % maximum number of Q-z points on one input line
    fprintf(F,'%s %s %s','*LAYER CODE:',element.model_axial{i},'for Q-z'); % Text to identify soil layer
    fprintf(F,'\n');
    if strcmp(plug_unplug.comp_index,'Unplugged')
        Q_rel = 1;
    elseif strcmp(plug_unplug.comp_index,'Plugged')
        Q_rel = pile.diameter^2/(pile.diameter^2-(pile.diameter-2*element.thickness(end))^2);
    end
    if Q_rel>=10
		fprintf(F,'%s %2.f  %5.2f %5.2f %6.2f','SOIL         SLOC   ',size(zQ.top,2),abs(element.level(end-1,1)-pile.head),abs(element.level(end-1,2)-pile.head),Q_rel); % Text to inform SACS that the following line contains user defined Q-z data. 
	else
		fprintf(F,'%s %2.f  %5.2f %5.2f %5.2f','SOIL         SLOC   ',size(zQ.top,2),abs(element.level(end-1,1)-pile.head),abs(element.level(end-1,2)-pile.head),Q_rel); % Text to inform SACS that the following line contains user defined Q-z data. 
	end
    if size(zQ.top,2) > 30 || size(zQ.bottom,2) > 30
        error('For use of more than 30 points on the t-z curve, you must specify the number in the header')
    end
    fprintf(F,'\n');
    for ii = 1:size(Qzdata,1)
            for jj = 1:size(Qzdata,2)
                if floor(jj/2) == jj/2 % even number, is valid for all z-values (lower precision is needed here compared to t-values)
                    if Qzdata(ii,jj) >= 10 % number is bigger than 10 (meaning there is room for one less digit after the comma)
                        indicator(ii,jj) = 1; 
                    else % number is smaller than 10
                        indicator(ii,jj) = 2;
                    end
                else % uneven number, is valid for all t-values (see comments further up)
                    if Qzdata(ii,jj) >= 10
                        indicator(ii,jj) = 3; 
                    else
                        indicator(ii,jj) = 4;
                    end
                end
            end
    end
    form = {'%4.1f ';'%4.2f ';'%6.3f ';'%6.4f '}; % form indicator is connected to numbers 1-4 shown above
    for nsec = 1:ceil(size(Qzdata,2)/(2*n_max)); % number of sections to print, maximum n_max t-z pairs in each section
        fprintf(F,'SOIL         T-Z ');
        form_indicator = [form{indicator(i,1+2*n_max*(nsec-1):min(2*n_max+2*n_max*(nsec-1),size(Qzdata,2)))}]; % chooses format for each line of t-z value pairs
        fprintf(F,form_indicator, Qzdata(i,1+2*n_max*(nsec-1):min(2*n_max+2*n_max*(nsec-1),size(Qzdata,2)))); % prints the t-z value pairs
        fprintf(F,'\n');
    end
end
if settings.torsion
    %--------------------------------------------------------------------------
    %% Torsion
    %--------------------------------------------------------------------------
    fprintf(F,'********************************************************************************\n');
    fprintf(F,'*TORSION');
    fprintf(F,'\n********************************************************************************\n');
    fprintf(F,'%s             %10.1f%s','SOIL TORSION HEAD',data.torsion,'TYP1');
	fprintf(F,'\n');
end
end
if settings.lateral_loading
%--------------------------------------------------------------------------
%% p-y curves
%--------------------------------------------------------------------------
% FORMAT
% Putting p.top, p.bottom and y into two matrices, py.top and py.bottom,
% where the format is [p1.1 y1.1  p1.2 y1.2  p1.3 y1.3  p1.4 y1.4  p1.5 y1.5 .. ;
%                      p2.1 y2.1  p2.2 y2.2  p2.3 y2.3  p2.4 y2,4  p2.5 y2.5 .. ; ..]
% ".." resembles the possibility that there can be more than more than five
% points per depth interval and more than two depth intervals
i = 1;
m = 1;
while i<=size(p.top,1)
    for j=1:size(y.top,2);
        % Top nodes
        py.top(m,j*2-1)     = p.top(i,j)/1E2; % to go from [kN/m] to [kN/cm]
        py.top(m,j*2)       = y.top(i,j)*1E2; % to go from [m] to [cm]
        % Bottom nodes
        py.bottom(m,j*2-1)  = p.bottom(i,j)/1E2; % to go from [kN/m] to [kN/cm]
        py.bottom(m,j*2)    = y.bottom(i,j)*1E2; % to go from [m] to [cm]
    end
    m = m + 1;
    i = i + 1;
end
% PRINT
fprintf(F,'********************************************************************************\n');
fprintf(F,'*P-Y CURVES');
fprintf(F,'\n********************************************************************************\n');
fprintf(F,'%s%3.0f        %5.1f       %s','SOIL LATERAL HEAD',2*length(element.level),pile.diameter*100,'TYP1'); % p-y curves are defined on top and bottom node of each element
fprintf(F,'\n');
for i = 1:size(element.level,1)-1                         %loop to print p-y curves
    for top_bottom = 1:2
        if top_bottom == 1
            pydata = py.top;
        else
            pydata = py.bottom;
        end
        n_max = 5; % maximum number of p-y points on one input line
        form = {'%5.2f ';'%5.1f '}; % format for p-y points
        indicator = pydata >= 100; % Boolean indicator to determine what format to use
        fprintf(F,'%s %s %s','*LAYER CODE:',element.model_py{i},'for lateral resistance and p-y'); % Text to identify soil layer
        fprintf(F,'\n');
        fprintf(F,'%s %2.f','SOIL         SLOCSM ',size(y.top,2)); % Text to inform SACS that the following line contains user defined p-y data
        if size(y.top,2) > 30 || size(y.bottom,2) > 30
            error('For use of more than 30 points on the p-y curve, you must specify the number in the header')
        end
        if top_bottom == 1
            % the print of the depth interval for each set of p-y values
            % can cause problems if depth is greater than 100 m - this is
            % because SACS only reads from certain columns and three
            % digits before the comma, e.g. 100.50, will not fit inside the
            % before mentioned columns
            if lin_con == 2 % curves are constant in the specified interval
                fprintf(F,'  %5.2f %5.2f', abs(element.level(i,1)-pile.head),abs(element.level(i,1)+(element.level(i,2)-element.level(i,1))/2-pile.head)); % Upper depth limit for p-y curve
            else % if not constant, then always varying linearly
                fprintf(F,'  %5.2f', abs(element.level(i,1)-pile.head)); % Upper depth limit for p-y curve
            end
        else
            if lin_con == 2 % curves are constant in the specified interval
                fprintf(F,'  %5.2f %5.2f', abs(element.level(i,1)+(element.level(i,2)-element.level(i,1))/2-pile.head),abs(element.level(i,2)-pile.head)); % Lower depth limit for p-y curves
            else % if not constant, then always varying linearly
                fprintf(F,'  %5.2f', abs(element.level(i,2)-pile.head)); % Upper depth limit for p-y curve
            end
        end
        fprintf(F,'\n');
        for nsec = 1:ceil(size(pydata,2)/(2*n_max)); % number of sections to print, maximum n_max p-y pairs in each section
            fprintf(F,'SOIL         P-Y  ');
            form_indicator = [form{indicator(i,1+2*n_max*(nsec-1):min(2*n_max+2*n_max*(nsec-1),size(pydata,2)))+1}]; % chooses format for each line of p-y value pairs
            fprintf(F,form_indicator, pydata(i,1+2*n_max*(nsec-1):min(2*n_max+2*n_max*(nsec-1),size(pydata,2)))); % prints the p-y value pairs
            fprintf(F,'\n');
        end
    end
end
end
fprintf(F,'********************************************************************************\n');
fprintf(F,'END');
fclose(F);