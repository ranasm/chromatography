%@param filepath Full filepath of ParentFraction.txt file
%@param rootpath Base folder filepath to output screenshot and text files

function UncorrectedBldHillFitPlotWrite(uncorrectedbloodpath,parentfractionpath)

    [rootpath,filename,~] = fileparts(uncorrectedbloodpath); %Store only the folder name, not full path

    uncorrectedblood=dlmread(uncorrectedbloodpath,'',1,0);
    
    uncorrectedblood(:,1)=uncorrectedblood(:,1)/60;
    time=uncorrectedblood(:,1);   %sampling times
    
    pfmetabolites=dlmread(parentfractionpath,'',1,0)
    pfmetabolites(:,2)=pfmetabolites(:,2)./100;

    'Fitting metabolites with Hill Fun'
    met0=[0.08;4.9;21.2];
    metlb=[0;0;0];            %mx0-4*stddev
    metub=[2;10;250];
    oo = optimset('MaxFunEvals', 5000,'MaxIter',3500);
    [metp,res,iterations]=lsqcurvefit(@hill_fun,met0,pfmetabolites(:,1),pfmetabolites(:,2),metlb,metub,oo);
    metabolite_fit=hill_fun(metp,time);
    plasmacor=uncorrectedblood(:,2).*metabolite_fit;
    
    %Plot
    figure('Name',sprintf('%s - Hill Function ',filename),'NumberTitle','off')
    plot(time,uncorrectedblood(:,2),'color','red','LineWidth',1.5);
    hold on;
    plot(time,plasmacor,'color','blue','LineWidth',1.5)
    legend({'Non Corrected Blood Plasma',' Corrected Blood Plasma'},'Location','Best')
    xlabel('Time (minutes)')
    ylabel('Radioactivity')
    title('Plasma Total Metabolite Corrected - Hill Function');
    
    %Screenshot
    set(gca,'fontsize', 16);
    saveas(gca,[rootpath sprintf('/%s_PlasmaTotalMetaboliteCorrected-Hill.png',filename)]);

    %Add parameter fits
    var=['a','b','c'];
    plen=length(metp);
    param = fullfile(rootpath,sprintf('/%s_PlasmaTotalMetaboliteCorrected-Hill_Fitted_Parameters.txt',filename));
    file = fopen(param,'w');
    fprintf(file,'Hill Fitted Parameters: \n');
    for j=1:plen
        fprintf(file,'%s = %f \n',var(j),metp(j))
    end
    fclose(file);
    
    %Write out to text file
    out = fullfile(rootpath,sprintf('/%s_PlasmaTotalMetaboliteCorrected-Hill.txt',filename));
    fileID = fopen(out,'w');
    fprintf(fileID,'sample-time[seconds]          whole-blood[kBq/cc]\n');
    Len=length(time);
    for i=1:Len
        fprintf(fileID,'%f,          %f \n',time(i)*60,plasmacor(i)); %Convert Back into seconds
    end
    fclose(fileID);
    
end

