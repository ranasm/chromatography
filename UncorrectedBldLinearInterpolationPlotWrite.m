%@param filepath Full filepath of ParentFraction.txt file
%@param rootpath Base folder filepath to output screenshot and text files

function UncorrectedBldLinearInterpolationPlotWrite(uncorrectedbloodpath,parentfractionpath)

    [rootpath,filename,~] = fileparts(uncorrectedbloodpath); %Store only the folder name, not full path

    uncorrectedblood=dlmread(uncorrectedbloodpath,'',1,0);
    
    uncorrectedblood(:,1)=uncorrectedblood(:,1)/60;
    time=uncorrectedblood(:,1);   %sampling times
    
    pfmetabolites=dlmread(parentfractionpath,'',1,0)
    pfmetabolites(:,2)=pfmetabolites(:,2)./100;
    
    %Fit Values
    interpolate=interp1(pfmetabolites(:,1),pfmetabolites(:,2),time,'linear','extrap'); 
    plasmacor=uncorrectedblood(:,2).*interpolate;
    
    %Plot 
    figure('Name',sprintf('%s - Interpolated',filename),'NumberTitle','off')
    plot(time,uncorrectedblood(:,2),'color','red','LineWidth',1.5);
    hold on;
    plot(time,plasmacor,'color','blue','LineWidth',1.5)
    legend({'Non Corrected Blood Plasma',' Corrected Blood Plasma'},'Location','Best')
    xlabel('Time (minutes)')
    ylabel('Radioactivity')
    title('Plasma Total Metabolite Corrected - Interpolated');
    
    set(gca,'fontsize', 16);

    %Screenshots
    saveas(gca,[rootpath sprintf('/%s_PlasmaTotalMetaboliteCorrected-Interpolated.png',filename)]);
    
    %Write to textfile
    out = fullfile(rootpath,sprintf('/%s_PlasmaTotalMetaboliteCorrected-Interpolated.txt',filename));
    fileID = fopen(out,'w');
    fprintf(fileID,'sample-time[seconds]          whole-blood[kBq/cc]\n');
    Len=length(time);
    for i=1:Len
        fprintf(fileID,'%f,          %f \n',time(i)*60,plasmacor(i));
    end
    fclose(fileID);


end

