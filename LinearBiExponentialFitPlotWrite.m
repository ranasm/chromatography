%@param filepath Full filepath of ParentFraction.txt file
%@param rootpath Base folder filepath to output screenshot and text files

function LinearBiExponentialFitPlotWrite(filepath,rootpath)

    metabolites=dlmread(filepath,'',1,0);
    time=0:1:metabolites(end,1);
    
    metabolites(:,2)=metabolites(:,2)./100;

    %Linear + Bi-Exponential
    var=["A1",'A3','lamda1','lamda2','tau'];
    tau=metabolites(2,1);
    met0=[-0.001; 0; 0.007367; 0.0003396;tau];
    metlb=[-1; 0; 0; 0;tau];
    metub=[0; 1; 1; 1;tau];
    oo = optimset('MaxFunEvals', 15000,'MaxIter',13500);
    [metp,res,iterations]=lsqcurvefit(@linear_biexp_fun,met0,metabolites(:,1),metabolites(:,2),metlb,metub,oo);
    metabolite_fit=linear_biexp_fun(metp,time);
    
    %Plot
    figure('Name','Fitted Linear + Bi-Exponential ','NumberTitle','off')
    plot(metabolites(:,1),metabolites(:,2),'*','color','red');
    hold on;
    plot(time,metabolite_fit,'color','blue');
    ylim([0 1])
    legend('samples','fitted Linear+Biexp')
    title('Fitted Linear + Bi-Exponential');
    set(gca,'fontsize', 16);
    [~,name,~] = fileparts(rootpath)  % This gives name of folder that contains all .D files

    %Screenshot
    saveas(gca,[rootpath sprintf('/proc_%s_linearbiexp_fitted_ParentFraction.png',name)]);
    biexpres=[time;metabolite_fit]';

    %Add parameter fits
    plen=length(metp);
    param = fullfile(rootpath,sprintf('/proc_%s_linearbiexp_fitted_ParentFraction_FitParameters.txt',name));
    file = fopen(param,'w');
    fprintf(file,'Linear+Biexp Fitted Parameters: \n');
    fprintf(file,'A2 = 1 \n');
    for j=1:plen
        fprintf(file,'%s = %f \n',var(j),metp(j));
    end
    fclose(file);

    %Write out text file
    out = fullfile(rootpath,sprintf('/proc_%s_linearbiexp_fitted_ParentFraction.txt',name))
    fileID = fopen(out,'w');
    fprintf(fileID,'Time              Percent Fraction \n');
    Len=length(biexpres);
    for i=1:Len
        fprintf(fileID,'%f,          %f \n',biexpres(i,1)*60,biexpres(i,2));
    end
    fclose(fileID);

end

