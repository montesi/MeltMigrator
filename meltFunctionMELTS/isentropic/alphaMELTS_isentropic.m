%% alphaMELTS Isentropic Mantle Upwelling & Melting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This script runs batch calculation for isentropic mantle upwelling & 
% melting using alphaMELTS with a sets of mantle potential temperature,
% and compiles the output data into a look-up table for later 
% interpolation.
%
% To have a control on the potential temperature, we start the calculation 
% at the surface with the desired potential temperature, and increases
% pressure to ~45000 [bar] (the maximum pressure alphaMELTS can handle), 
% with liquid phase supressed. This gives us initial conditions and 
% composition for mantle upwelling. Starting from there, we turn on the 
% liquid phase and decrease pressure adiabatically and track the evolution
% of the system. 
%
% Hailogn Bai & Laurent Montesi
% 2015-05-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization

clear;
T0=[1450:-5:1100]; % Mantle potential temperature
IsentropicDirectory='E:\Scripts\MeltMigrator\meltFunctionMELTS\isentropic';
cd(sprintf('%s',IsentropicDirectory));

Switch_BatchFilesCreated=1;
Switch_alphaMELTsRunDone=1;

%% Create batch file for alphaMELTS calculation for each potential temperature

if ~Switch_BatchFilesCreated
    cd batch; % Batch file directory
    for ind=1:numel(T0);
        fid0=fopen('batch.txt','r+'); % change second line of batch.txt
        count=0;
        strg=[];
        while ~feof(fid0);
            tline=fgetl(fid0);
            count=count+1;
            strg{count}=tline;
            if count==4;
                strg{count}=strrep(strg{count},'1375',sprintf('%g',T0(ind)));
            end
        end
        fclose(fid0);
        fid=fopen(sprintf('batch_%g.txt',T0(ind)),'wt');
        for j=1:count;
            fprintf(fid,'%s\n',strg{j});
        end
        fclose(fid);
    end
end

%% Batch alphaMELTS

if ~Switch_alphaMELTsRunDone
    SuccessFlag_All=[];
    tic;
    for i=[1:numel(T0)];  
        cd(sprintf('%s',IsentropicDirectory));
        % Batch alphaMELTS calculation
        mkdir(sprintf('output\\%g',T0(i)));
        cd(sprintf('%s\\output\\%g',IsentropicDirectory,T0(i))); % Output directory for each calculation
        [status readout]=system(sprintf('run_alphamelts.command -f %s\\isentropic_melt_env.txt -p %s\\output\\%g -b %s\\batch\\batch_%g.txt',IsentropicDirectory,IsentropicDirectory,T0(i),IsentropicDirectory,T0(i)),'-echo');
        fileID=fopen('readout_file.txt','w');
        fprintf(fileID,'%s',readout);
        fclose(fileID);

        % Determine if the calculation is successful
        fmelts=fopen('melts_file.melts','r+');
        strind=[];
        pressureline=['This is the line where the initial pressure is supposed to be given'];
        while (isempty(strind))&&(~feof(fmelts));
            tline=fgetl(fmelts);
            strind=strfind(tline,'Initial Pressure');
            if ~isempty(strind);
                pressureline=tline;
            end
        end
        fclose(fmelts);
        if (pressureline(19)=='1')&&(pressureline(20)=='.'); % If reach the pressure of 1.00 bar, the calculation is finished
            flag_finish=1;
        else
            flag_finish=0;
        end
        SuccessFlag_All=[SuccessFlag_All;flag_finish];

        % If error did occur, use the saved .melts file to continue the calculation
        % (Usually, the second calculation will make it successfully to the surface, if not, a third calculation may be needed)
        if flag_finish==0;
            cd(sprintf('%s\\batch',IsentropicDirectory));
            fid0=fopen('batch_error.txt','r+'); % change second line of batch_error.txt
            count=0;
            strg=[];
            while ~feof(fid0);
                tline=fgetl(fid0);
                count=count+1;
                strg{count}=tline;
                if count==2;
                    strg{count}=strrep(strg{count},'1375',sprintf('%g',T0(i)));
                end
            end
            fclose(fid0);
            fid=fopen(sprintf('batch_error_%g.txt',T0(i)),'wt');
            for j=1:count;
                fprintf(fid,'%s\n',strg{j});
            end
            fclose(fid);
            [status readout]=system(sprintf('run_alphamelts.command -f %s\\isentropic_melt_env_error.txt -p %s\\output\\%g\\continued -b %s\\batch\\batch_error_%g.txt',IsentropicDirectory,IsentropicDirectory,T0(i),IsentropicDirectory,T0(i)),'-echo');
            cd(sprintf('%s\\output\\%g',IsentropicDirectory,T0(i)));
            fileID=fopen('readout_file_continued.txt','w');
            fprintf(fileID,'%s',readout);
            fclose(fileID);
        end
    end

    cd(sprintf('%s',IsentropicDirectory));
    save('SuccessFlag_All');
    t=toc;
    display(sprintf('Processing time: %s',secs2hms(t))); % Display the time it took to finish the batch calculation
end

%% Read data from alphaMELTS output

cd(sprintf('%s',IsentropicDirectory));
load('SuccessFlag_All.mat'); 
M_total=100.139;
T_ALL=[];
M_ALL=[];
F_ALL=[];

for i=[1:numel(T0)];
    cd(sprintf('%s\\output\\%g',IsentropicDirectory,T0(i)));
    
    % Choose the output file of interest
    fid=fopen('System_main_tbl.txt');
    C=textscan(fid,'%f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f','HeaderLines',4,'Delimiter',' ');
    % (1)Pressure (2)Temperature (3)mass (4)F (5)phi (6)S (7)H (8)V (9)Cp (10)dVdP*10^6 (11)dVdT*10^6 (12)fO2(absolute) (13)fO2(absolute) (14)rhol (15)rhos (16)viscosity (17)aH2O (18)chisqr
    
%     fid=fopen('Solid_comp_tbl.txt');
%     C=textscan(fid,'%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f','HeaderLines',4,'Delimiter',' ');
%     % (1)Pressure (2)Temperature (3)mass (4)SiO2 (5)TiO2 (6)Al2O3 (7)Fe2O3 (8)Cr2O3 (9)FeO (10)MnO (11)MgO (12)NiO (13)CaO (14)Na2O

%     fid=fopen('Phase_mass_tbl.txt');
%     C=textscan(fid,'%f %f %f %f %*f %*f %*f %*f %*f %*f','HeaderLines',4,'Delimiter',' ');
%     % (1)Pressure (2)Temperature (3)mass (4)liquid (5)olivine (6)garnet (7)orthopyroxene (8)clinopyroxene (9)feldspar (10)spinel

%     fid=fopen('Bulk_comp_tbl.txt');
%     C=textscan(fid,'%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f','HeaderLines',4,'Delimiter',' ');
%     % (1)Pressure (2)Temperature (3)mass (4)SiO2 (5)TiO2 (6)Al2O3(7)Fe2O3 (8)Cr2O3 (9)FeO (10)MnO (11)MgO (12)NiO (13)CaO (14)Na2O
    
%     Phase_main_tbl.txt, Liquid_comp_tbl.txt have irragule size and are difficult to read

    fclose(fid);
    
    if SuccessFlag_All(i)==0;
        cd(sprintf('%s\\output\\%g\\continued',IsentropicDirectory,T0(i)));
        fid=fopen('System_main_tbl.txt');
        D=textscan(fid,'%f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f','HeaderLines',4,'Delimiter',' ');
        fclose(fid);
        Palpha=[C{1};D{1}];
        Talpha=[C{2};D{2}];
        Malpha=[C{3};D{3}];
        Falpha=[C{4};D{4}];
    else
        Palpha=C{1};
        Talpha=C{2};
        Malpha=C{3};
        Falpha=C{4};
    end
    
    T_ALL=[T_ALL,Talpha];
    M_ALL=[M_ALL,Malpha];
    F_ALL=[F_ALL,Falpha];
end

FiALL=(M_total-M_ALL)./M_total;

%% Plot the final look-up table

[Tg,Pg]=meshgrid(T0,Palpha);

figure(1); clf;
pcolor(Tg,Pg,T_ALL);shading interp; colorbar;
% pcolor(Tg,Pg,TALL-Tg);shading interp; colorbar;
set(gca,'ydir','reverse');
xlabel('Mantle Potential Temperature [\circC]');
ylabel('Pressure [bar]');
title('Temperature [\circC]');
set(gca,'fontsize',14);
figure(2); clf;
pcolor(Tg,Pg,FiALL);shading interp; colorbar;
set(gca,'ydir','reverse');
xlabel('Mantle Potential Temperature [\circC]');
ylabel('Pressure [bar]');
title('Melt Fraction');
set(gca,'fontsize',14);
% set(gca,'clim',[0,0.25]);
