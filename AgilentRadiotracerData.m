classdef AgilentRadiotracerData < handle
    % Wrapper for chromatography class
    % Detailed explanation goes here
    
    properties
        initializetoolbox;
        radiotracerdata;
    end
    
    methods
        function obj = AgilentRadiotracerData()
            %   Constructor
            %   This will prompt user to load in select agilent folder
            obj.initializetoolbox=Chromatography();
            obj.radiotracerdata=obj.initializetoolbox.import('.D');
        end
        
        function DecayCorrectRadiotracer(obj,typeofisotope)
            if typeofisotope=='F-18'
                obj.DecayCorrectF18();
                
            elseif typeofisotope=='C-11'
                obj.DecayCorrectC11();
            else 
                error('Could not apply proper decay correction to specified radiotracer.');
            end
        end
        
         %-------------------------
         % Refer to http://www.turkupetcentre.net/petanalysis/decay.html
         % for derivation of correction function to apply 
         %-------------------------
        function DecayCorrectC11(obj)
            lamda=.034;
            lamdatimeproduct=lamda.*(obj.radiotracerdata.time);
            exponentterm=exp(lamdatimeproduct);
            decaycorrectedvalues=obj.radiotracerdata.tic.values.*exponentterm;
            obj.radiotracerdata.tic.values=decaycorrectedvalues;
        end
        
        function DecayCorrectF18(obj)
            lamda=.006;
            lamdatimeproduct=lamda.*(obj.radiotracerdata.time);
            exponentterm=exp(lamdatimeproduct);
            decaycorrectedvalues=obj.radiotracerdata.tic.values.*exponentterm;
            obj.radiotracerdata.tic.values=decaycorrectedvalues;
        end
        
        function ApplySmoothingandAsymmetry(obj,smooth,asymmetry)
            obj.radiotracerdata=obj.initializetoolbox.smooth(obj.radiotracerdata,...
                'samples','all',...
                'ions','tic',...
                'smoothness',smooth,...
                'asymmetry',asymmetry)
        end
        
        function ShiftTracerDataAboveBaseline(obj)
            translationamount=abs(min(obj.radiotracerdata.tic.values));
            obj.radiotracerdata.tic.values=obj.radiotracerdata.tic.values+translationamount;
        end
        
        %% To do
        function WriteTracerDatatoText(obj)
            disp('');
        end
     
    end
end 