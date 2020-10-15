%@param filepath Full filepath of ParentFraction.txt file
%@param rootpath Base folder filepath to output screenshot and text files

function HillFitPlotWrite(filepath,rootpath)

    metabolites=dlmread(filepath,'',1,0);
    time=0:1/60:metabolites(end,1);
    
    metabolites(:,2)=metabolites(:,2)./100;
    
    %Hill Function
    'Fitting metabolites with Hill Fun'
    met0=[0.08;4.9;21.2];
    metlb=[0;0;0];            %mx0-4*stddev
    metub=[2;10;250];
    oo = optimset('MaxFunEvals', 5000,'MaxIter',3500);
    [metp,res,iterations]=lsqcurvefit(@hill_fun,met0,metabolites(:,1),metabolites(:,2),metlb,metub,oo);
    metabolite_fit=hill_fun(metp,time);
    
    %Plot
    figure('Name','Fitted Hill Function ','NumberTitle','off')
    plot(metabolites(:,1),metabolites(:,2),'*','color','red');
    hold on;
    plot(time,metabolite_fit,'color','blue')
    ylim([0 1])
    title('Fitted Hill Function');
    legend('samples','fitted Hill')
    set(gca,'fontsize', 16);
    [~,name,~] = fileparts(rootpath)  % This gives name of folder that contains all .D files
    
    %Screenshot
    saveas(gca,[rootpath sprintf('/proc_%s_hill_fitted_ParentFraction.png',name)]);
    hillfitres=[time;metabolite_fit]';

    %Write fitted values to text file
    out = fullfile(rootpath,sprintf('/proc_%s_Hill_Fitted_ParentFraction.txt',name));
    fileID = fopen(out,'w');
    fprintf(fileID,'Time              Percent Fraction \n');
    Len=length(hillfitres);
    for i=1:Len
        fprintf(fileID,'%f,          %f \n',hillfitres(i,1)*60,hillfitres(i,2));
    end
    fclose(fileID);    
    
    %Write parameters to text file
    var=['a','b','c'];
    plen=length(metp);
    param = fullfile(rootpath,sprintf('/proc_%s_Hill_Fitted_ParentFraction_FitParameters.txt',name));
    file = fopen(param,'w');
    fprintf(file,'Hill Fitted Parameters: \n');
    for j=1:plen
        fprintf(file,'%s = %f \n',var(j),metp(j));
    end
    fclose(file);

end

