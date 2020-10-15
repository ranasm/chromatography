%@param filepath Full filepath of ParentFraction.txt file
%@param rootpath Base folder filepath to output screenshot and text files

function UncorrectedBldLinearBiExponentialFitPlotWrite(uncorrectedbloodpath,parentfractionpath)
    
    [rootpath,filename,~] = fileparts(uncorrectedbloodpath); %Store only the folder name, not full path
    
    uncorrectedblood=dlmread(uncorrectedbloodpath,'',1,0);
    
    uncorrectedblood(:,1)=uncorrectedblood(:,1)/60;
    time=uncorrectedblood(:,1);   %sampling times
    
    pfmetabolites=dlmread(parentfractionpath,'',1,0)
    pfmetabolites(:,2)=pfmetabolites(:,2)./100;
    
    % Linear + Biexponential Values
    tau=pfmetabolites(2,1);
    met0=[-0.001; 0; 0.007367; 0.0003396;tau];
    metlb=[-1; 0; 0; 0;tau];
    metub=[0; 1; 1; 1;tau];
    oo = optimset('MaxFunEvals', 15000,'MaxIter',13500);
    [metp,res,iterations]=lsqcurvefit(@linear_biexp_fun,met0,pfmetabolites(:,1),pfmetabolites(:,2),metlb,metub,oo);
    metabolite_fit=linear_biexp_fun(metp,time);
    plasmacor=uncorrectedblood(:,2).*metabolite_fit;
    
    %Plot
    figure('Name',sprintf('%s - Fitted Linear + BiExponential ',filename),'NumberTitle','off')
    plot(time,uncorrectedblood(:,2),'color','red','LineWidth',1.5);
    hold on;
    plot(time,plasmacor,'color','blue','LineWidth',1.5)
    legend({'Non Corrected Blood Plasma',' Corrected Blood Plasma'},'Location','Best')
    title('Plasma Total Metabolite Corrected - Fitted Linear + BiExponential');
    set(gca,'fontsize', 16);
    xlabel('Time (minutes)')
    ylabel('Radioactivity')
    
    %Screenshot
    saveas(gca,[rootpath sprintf('/%s_PlasmaTotalMetaboliteCorrected-Linear+BiExponential.png',filename)]);
    
    %Write out fitted values
    plen=length(metp);
    var=["A1",'A3','lamda1','lamda2','tau'];
    param = fullfile(rootpath,sprintf('/%s_PlasmaTotalMetaboliteCorrected-Linear+BiExponential_Fitted_Parameters.txt',filename));
    file = fopen(param,'w');
    fprintf(file,'Linear+BiExponential Fitted Parameters: \n');
    fprintf(file,'A2 = 1 \n');
    for j=1:plen
        fprintf(file,'%s = %f \n',var(j),metp(j));
    end
    fclose(file);

    %Write out corrected data
    out = fullfile(rootpath,sprintf('/%s_PlasmaTotalMetaboliteCorrected-Linear+BiExponential.txt',filename));
    fileID = fopen(out,'w');
    fprintf(fileID,'sample-time[seconds]          whole-blood[kBq/cc]\n');
    Len=length(time);
    for i=1:Len
        fprintf(fileID,'%f,          %f \n',time(i)*60,plasmacor(i));
    end
    fclose(fileID);

end

