function MeltFraction=meltFunctionMELTS(z,Temperature_COMSOL)

%% Generate the look-up table

% % use alphaMELTS data to modify thermal structure
PreviousDirectory=pwd;
cd F:\Migration\alphaMELTS\isentropic\;
Temperature_Potential=[1450:-5:1100]; % potential temperature [degC]
Mass_Total=100.139;
Density_Mantle=3300;
g=9.8;

Flag_alphaMELTS=load('flag_all.mat'); 
Flag_alphaMELTS=Flag_alphaMELTS.flag_all;
Temperature_alphaMELTS=[];
Mass_alphaMELTS=[];

for i_PoT=[1:numel(Temperature_Potential)];
    cd(sprintf('F:\\Migration\\alphaMELTS\\isentropic\\output\\%g',Temperature_Potential(i_PoT)));
    
    % choose the output file of interest
    Output_alphaMELTS=fopen('System_main_tbl.txt');
    Data_alphaMELTS=textscan(Output_alphaMELTS,'%f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f','HeaderLines',4,'Delimiter',' ');
    % (1)Pressure (2)Temperature (3)mass (4)F (5)phi (6)S (7)H (8)V (9)Cp (10)dVdP*10^6 (11)dVdT*10^6 
    % (12)fO2(absolute) (13)fO2(absolute) (14)rhol (15)rhos (16)viscosity (17)aH2O (18)chisqr
    fclose(Output_alphaMELTS);
    
    if Flag_alphaMELTS(i_PoT)==0;
        cd(sprintf('F:\\Migration\\alphaMELTS\\isentropic\\output\\%g\\continued',Temperature_Potential(i_PoT)));
        Output_alphaMELTS=fopen('System_main_tbl.txt');
        Data_alphaMELTS_Continued=textscan(Output_alphaMELTS,'%f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f','HeaderLines',4,'Delimiter',' ');
        fclose(Output_alphaMELTS);
        Pressure_alphaMELTS=[Data_alphaMELTS{1};Data_alphaMELTS_Continued{1}]./1e4;
        Temperature_alphaMELTS_Column=[Data_alphaMELTS{2};Data_alphaMELTS_Continued{2}];
        Mass_alphaMELTS_Column=[Data_alphaMELTS{3};Data_alphaMELTS_Continued{3}];
    else
        Pressure_alphaMELTS=Data_alphaMELTS{1}./1e4;
        Temperature_alphaMELTS_Column=Data_alphaMELTS{2};
        Mass_alphaMELTS_Column=Data_alphaMELTS{3};
    end    
    Temperature_alphaMELTS=[Temperature_alphaMELTS,Temperature_alphaMELTS_Column];
    Mass_alphaMELTS=[Mass_alphaMELTS,Mass_alphaMELTS_Column];
end

MeltFraction_alphaMELTS=(Mass_Total-Mass_alphaMELTS)./Mass_Total;
[GridTemperature,GridPressure]=meshgrid(Temperature_Potential,Pressure_alphaMELTS);

%% Use look-up table to modify COMSOL temperature

Temperature_Modified=zeros(size(Temperature_COMSOL));
MeltFraction_Modified=zeros(size(Temperature_COMSOL));
z2P=@(z) Density_Mantle*g*abs(z)/1e6;

Pressure_COMSOL=z2P(z);

for ind_z=1:size(Temperature_COMSOL,3);
    for ind_y=1:size(Temperature_COMSOL(:,:,ind_z),2);
        Temperature_Temp=interp2(GridTemperature,GridPressure,Temperature_alphaMELTS,Temperature_COMSOL(:,ind_y,ind_z),Pressure_COMSOL(:,ind_y,ind_z));
        MeltFraction_Temp=interp2(GridTemperature,GridPressure,MeltFraction_alphaMELTS,Temperature_COMSOL(:,ind_y,ind_z),Pressure_COMSOL(:,ind_y,ind_z));
        Temperature_Modified(:,ind_y,ind_z)=Temperature_Temp;
        MeltFraction_Modified(:,ind_y,ind_z)=MeltFraction_Temp;
    end
end

indNaN=find(isnan(Temperature_Modified));
Temperature_Modified(indNaN)=Temperature_COMSOL(indNaN);
MeltFraction_Modified(isnan(MeltFraction_Modified))=0;
MeltFraction=MeltFraction_Modified;

cd(sprintf('%s',PreviousDirectory));

return