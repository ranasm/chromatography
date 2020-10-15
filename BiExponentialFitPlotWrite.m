%@param filepath Full filepath of ParentFraction.txt file
%@param rootpath Base folder filepath to output screenshot and text files

function BiExponentialFitPlotWrite(filepath,rootpath)

    metabolites=dlmread(filepath,'',1,0);
    time=0:1:metabolites(end,1);
    
    metabolites(:,2)=metabolites(:,2)./100;

    %Bi Exponential
    var=['a','b','c','d'];
    met0=[1;0.01;0.15;0.6];
    metlb=[0;0;0;0];
    metub=[3;2;2;2];
    oo = optimset('MaxFunEvals', 5000,'MaxIter',3500);
    [metp,res,iterations]=lsqcurvefit(@biexp_fun,met0,metabolites(:,1),metabolites(:,2),metlb,metub,oo);
    metabolite_fit=biexp_fun(metp,time);
    
    %Plot
    figure('Name','Fitted Bi-Exponential ','NumberTitle','off')
    plot(metabolites(:,1),metabolites(:,2),'*','color','red');
    hold on;
    plot(time,metabolite_fit,'color','blue');
    ylim([0 1])
    legend('samples','fitted Biexp')
    title('Fitted Bi-Exponential');
    set(gca,'fontsize', 16);
    [~,name,~] = fileparts(rootpath)  % This gives name of folder that contains all .D files
    
    %Screenshot
    saveas(gca,[rootpath sprintf('/proc_%s_biexp_fitted_ParentFraction.png',name)]);
    biexpres=[time;metabolite_fit]';

    %Add parameter fits
    plen=length(metp);
    param = fullfile(rootpath,sprintf('/proc_%s_biexp_fitted_ParentFraction_FitParameters.txt',name));
    file = fopen(param,'w');
    fprintf(file,'Biexp Fitted Parameters: \n');
    for j=1:plen
        fprintf(file,'%s = %f \n',var(j),metp(j));
    end
    fclose(file);

    %Write out fitted values
    out = fullfile(rootpath,sprintf('/proc_%s_BiExp_ParentFraction.txt',name))
    fileID = fopen(out,'w');
    fprintf(fileID,'Time              Percent Fraction\n');
    Len=length(biexpres)
    for i=1:Len
        fprintf(fileID,'%f,          %f \n',biexpres(i,1)*60,biexpres(i,2));
    end
    fclose(fileID);
end

