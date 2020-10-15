%@param filepath Full filepath of ParentFraction.txt file
%@param rootpath Base folder filepath to output screenshot and text files

function UncorrectedBldBiExponentialFitPlotWrite(uncorrectedbloodpath,parentfractionpath)
    
    [rootpath,filename,~] = fileparts(uncorrectedbloodpath); %Store only the folder name, not full path

    uncorrectedblood=dlmread(uncorrectedbloodpath,'',1,0);
    
    uncorrectedblood(:,1)=uncorrectedblood(:,1)/60;
    time=uncorrectedblood(:,1);   %sampling times
    
    pfmetabolites=dlmread(parentfractionpath,'',1,0)
    pfmetabolites(:,2)=pfmetabolites(:,2)./100;
   
    %Bi-Exponential
    met0=[1;0.01;0.15;0.6];
    metlb=[0;0;0;0];
    metub=[3;2;2;2];
    oo = optimset('MaxFunEvals', 5000,'MaxIter',3500);
    [metp,res,iterations]=lsqcurvefit(@biexp_fun,met0,pfmetabolites(:,1),pfmetabolites(:,2),metlb,metub,oo);
    metabolite_fit=biexp_fun(metp,time);
    plasmacor=uncorrectedblood(:,2).*metabolite_fit;
        
    %Plot
    figure('Name',sprintf('%s - Bi-Exponential Fit ',filename),'NumberTitle','off')
    plot(time,uncorrectedblood(:,2),'color','red','LineWidth',1.5);
    hold on;
    plot(time,plasmacor,'color','blue','LineWidth',1.5)
    legend({'Non Corrected Blood Plasma',' Corrected Blood Plasma'},'Location','Best')
    title('Plasma Total Metabolite Corrected - Fitted BiExponential');
    set(gca,'fontsize', 16);
    xlabel('Time (minutes)')
    ylabel('Radioactivity')
    
    %Screenshot
    saveas(gca,[rootpath sprintf('/%s_PlasmaTotalMetaboliteCorrected-Bi-Exponential.png',filename)]);

    %Add parameter fits
    var=['a','b','c','d'];
    plen=length(metp);
    param = fullfile(rootpath,sprintf('/%s_PlasmaTotalMetaboliteCorrected-Bi-Exponential_Fitted_Parameters.txt',filename));
    file = fopen(param,'w');
    fprintf(file,'Biexp Fitted Parameters: \n');
    for j=1:plen
        fprintf(file,'%s = %f \n',var(j),metp(j));
    end
    fclose(file);

    out = fullfile(rootpath,sprintf('/%s_PlasmaTotalMetaboliteCorrected-Bi-Exponential.txt',filename));
    fileID = fopen(out,'w');
    fprintf(fileID,'sample-time[seconds]          whole-blood[kBq/cc]\n');
    Len=length(time)
    for i=1:Len
        fprintf(fileID,'%f,          %f \n',time(i)*60,plasmacor(i));
    end
    fclose(fileID);



end

