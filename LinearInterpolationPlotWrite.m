%@param filepath Full filepath of ParentFraction.txt file
%@param rootpath Base folder filepath to output screenshot and text files
function LinearInterpolationPlotWrite(filepath,rootpath)

    metabolites=dlmread(filepath,'',1,0);
    time=0:1/60:metabolites(end,1);

    %Interpolated Parent Function
    metabolites(:,2)=metabolites(:,2)./100;
    interpolate=interp1(metabolites(:,1),metabolites(:,2),time);
    
    %Plot
    figure('Name','Interpolated Parent Fraction ','NumberTitle','off')
    plot(time,interpolate);
    ylim([0 1])
    hold on;
    plot(metabolites(:,1),metabolites(:,2),'*');
    xlabel('Time (minutes)')
    ylabel('Parent Fraction')
    title('Interpolated Parent Fraction');
    set(gca,'fontsize', 16);
    [~,name,~] = fileparts(rootpath)  % This gives name of folder that contains all .D files
    
    %Screenshot
    saveas(gca,[rootpath, sprintf('/proc_%s_Interpolated_ParentFraction.png',name)]);
    interpres=[time;interpolate]';
    out = fullfile(rootpath,sprintf('/proc_%s_Interpolated_ParentFraction.txt', name));
    
    %Write out textfile
    fileID = fopen(out,'w');
    fprintf(fileID,'Time              Percent Fraction\n');
    Len=length(interpres);
    for i=1:Len
        fprintf(fileID,'%f,          %f \n',interpres(i,1)*60,interpres(i,2));
    end
    fclose(fileID);


end