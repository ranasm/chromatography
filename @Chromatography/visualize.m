% Method: visualize
%  -Plot chromatograms, baselines, and curve fitting results
%
% Syntax
%   fig = visualize(data)
%   fig = visualize(data, 'OptionName', optionvalue...)
%
% Options
%   'samples'   : 'all', [index]
%   'ions'      : 'all', 'tic', [index]
%   'baseline'  : 'on', 'off', 'corrected'
%   'peaks'     : 'on', 'off', 'residuals'
%   'layout'    : 'stacked', 'overlaid'
%   'scale'     : 'normalized', 'full'
%   'scaling'   : 'local', 'global'
%   'padding'   : value
%   'offset'    : value
%   'linewidth' : value
%   'xlim'      : 'auto', [xmin, xmax]
%   'ylim'      : 'auto', [ymin, ymax]
%   'legend'    : 'on', 'off'
%   'export'    : see MATLAB documentation on print functions
%   'colormap'  : 'parula', 'jet', 'hsv', 'hot', 'cool', 'spring', 'summer',
%                 'autumn', 'winter', 'gray', 'bone', 'copper', 'pink'
%
% Description
%   data        : data structure
%   'samples'   : row index of samples (default = 'all')
%   'ions'      : column index of ions (default = 'tic')
%   'baseline'  : display baseline results (default = 'off')
%   'peaks'     : display curve fitting results (default = 'off')
%   'layout'    : stacked or overlaid spectra (default = 'overlaid')
%   'scale'     : normalized or full scale (default = 'full')
%   'scaling'   : scaling relative to local or global maxima (default = 'global')
%   'padding'   : whitespace between axes and spectra (default = 0.05)
%   'offset'    : y-offset between spectra (default = 0.0);
%   'linewidth' : line width of spectra (default = 1.5)
%   'xlim'      : x-axis limits (default = 'auto')
%   'ylim'      : y-axis limits (default = 'auto')
%   'legend'    : display legend (default = 'off')
%   'export'    : options passed to the MATLAB print function (default = 'off')
%   'colormap'  : select colormap to use for plotting (default = 'parula')
%
% Examples
%   fig = obj.visualize(data, 'samples', [1:4], 'ions', 'tic')
%   fig = obj.visualize(data, 'layout', 'stacked', 'ions', 'all')
%   fig = obj.visualize(data, 'scale', 'normalized', 'legend', 'on')
%   fig = obj.visualize(data, 'samples', 1, 'xlim', [3, 27])
%   fig = obj.visualize(data, 'ions', 'all', 'colormap', 'winter', 'export', {'myfig', '-dtiff', '-r300'}

function varargout = visualize(obj, varargin)

% Check input
[data, options] = parse(varargin);

% Variables
samples = options.samples;
ions = options.ions;
scale =  options.scale;

% Initialize axes
options = plot_axes(obj, options, data);

% Initialize plot
for i = 1:length(samples)
    
    % Counter
    options.i = i;
    
    % Input values
    x = data(samples(i)).time;
    
    % Check ion options
    if isnumeric(ions)
        ions = 'xic';
    end
    
    switch ions
        
        case 'tic'
            
            % Input values
            y = data(samples(i)).tic.values;
                    
            % Check baseline
            if any(strcmpi(options.baseline, {'on', 'corrected'}))
                baseline = data(samples(i)).tic.baseline;
            else
                baseline = [];
            end
            
            % Check peaks
            if any(strcmpi(options.peaks, {'on', 'residuals'}))
                peaks = data(samples(i)).tic.peaks.fit;
            else
                peaks = [];
            end
            
            % Display names
            options.name(end+1) = {data(samples(i)).sample.name};
            options.name = (unique(options.name, 'stable'));
                
        case 'all'
            
            % Input values
            y = data(samples(i)).xic.values;
                
            % Check baseline
            if any(strcmpi(options.baseline, {'on', 'corrected'}))
                baseline = data(samples(i)).xic.baseline;
            else
                baseline = [];
            end
               
            % Check peaks
            if any(strcmpi(options.peaks, {'on', 'residuals'}))
                peaks = data(samples(i)).xic.peaks.fit;
            else
                peaks = [];
            end
            
            % Display names
            for j = 1:length(y(1,:))
                options.name{j} = strcat(num2str(data(samples(i)).mz(j)), ' m/z');
            end
                
            % Filter duplicate names
            options.name = (unique(options.name, 'stable'));
        
        case 'xic'
                
            % Input values
            y = data(samples(i)).xic.values(:, options.ions);
        
            % Check baseline
            if any(strcmpi(options.baseline, {'on', 'corrected'}))
                baseline = data(samples(i)).xic.baseline;
                
                if ~isempty(baseline)
                    baseline = baseline(:, options.ions);
                end
            else
                baseline = [];
            end
               
            % Check peaks
            if any(strcmpi(options.peaks, {'on', 'residuals'}))
                peaks = data(samples(i)).xic.peaks.fit;
                
                if ~isempty(peaks)
                    peaks = peaks(:, options.ions);
                end
            else
                peaks = [];
            end
        
            % Display names
            for j = 1:length(options.ions)
                options.name{j} = strcat(num2str(data(samples(i)).mz(options.ions(j))), ' m/z');
            end
        
            % Filter duplicate names
            options.name = (unique(options.name, 'stable'));
    end
    
    % Check baseline
    if ~isempty(baseline) && strcmpi(options.baseline, 'corrected');
        y = y - baseline;
    elseif isempty(baseline)
        options.baseline = 'off';
    end
    
    % Check peaks
    if isempty(peaks)
        options.peaks = 'off';
    end
    
    % Apply user input
    y = plot_scale(y, options);
    y = plot_layout(y, options);
    
    options = plot_xlim(x, options);
    options = plot_ylim(y, options);
    
    % Plot data
    switch version('-release')
        
        case '2014b'
            for j = 1:length(y(1,:))
                plot(x, y(:,j), ...
                    'parent', options.axes, ...
                    'linewidth', options.linewidth, ...
                    'displayname', options.name{j});
            end
            
        otherwise
            for j = 1:length(y(1,:))
                plot(x, y(:,j), ...
                    'parent', options.axes, ...
                    'linewidth', options.linewidth, ...
                    'linesmoothing', 'on',...
                    'displayname', options.name{j});
            end
    end
    
    % Check baseline options
    if strcmpi(options.baseline, 'on')
        
        if strcmpi(scale, 'normalized')
            baseline = (baseline - min(min(y))) / (max(max(y))-min(min(y)));
        end
        
        baseline = plot_layout(baseline, options);
        
        % Plot baseline
        plot(x, baseline, ...
            'parent', options.axes, ...
            'linewidth', options.linewidth, ...
            'color', [0.99,0.25,0.23]);
    end
end

% Display figure
set(options.figure, 'visible', 'on');

% Update plot properties
plot_update(options);

% Update axes position
set(options.empty, 'position', get(options.axes, 'position'));

% Export figure
if ~isempty(options.export)
    try
        disp('Rendering image, please wait...');
        print(options.figure, options.export{:});
        disp('Rendering complete!');
    catch
        disp('-Error reading print options, rendering image with default options...');
        print(options.figure, 'spectra', '-dpng', '-r300');
        disp('Rendering complete!');
    end
end

% Output
varargout{1} = options;

end


% Check scale options
function y = plot_scale(y, options)

% Determine scale
if strcmpi(options.scale, 'normalized') && strcmpi(options.scaling, 'local')
    y = bsxfun(@rdivide, bsxfun(@minus, y, min(y)), (max(y)-min(y)));

elseif strcmpi(options.scale, 'normalized') && strcmpi(options.scaling, 'global')
    y = bsxfun(@rdivide, bsxfun(@minus, y, min(min(y))), (max(max(y))-min(min(y))));

elseif strcmpi(options.scale, 'full') && strcmpi(options.layout, 'stacked')
    
    if strcmpi(options.scaling, 'local')
        y = bsxfun(@rdivide, bsxfun(@minus, y, min(min(y))), (max(max(y))-min(min(y))));
    end
end


end


% Check layout options
function y = plot_layout(y, options)

% Determine stacked layout
if strcmpi(options.layout, 'stacked') && length(options.samples) > 1
    
    % Determine scale
    if strcmpi(options.scale, 'normalized')
        
        % Calculate offset
        y = y - options.i * (1 + options.padding + options.offset);
    
    elseif strcmpi(options.scale, 'full')
        
        % Determine y-limits
        if isempty(options.ylimits)
            options.ylimits(1) = min(min(y));
            options.ylimits(2) = max(max(y));
        end
        
        % Calculate offset
        padding = options.ylimits(2) * options.padding;
        offset = options.ylimits(2) * options.offset;
        
        y = y - (options.i-1) * (options.ylimits(2) + padding + offset);
    end
    
elseif strcmpi(options.layout, 'overlaid') && length(options.samples) > 1
    
    % Determine scale
    if options.offset ~= 0 && strcmpi(options.scale, 'normalized')
        
        % Calculate offset
        y = y - options.i * options.offset;

    elseif options.offset ~= 0 && strcmpi(options.scale, 'full')

        % Determine y-limits
        if isempty(options.ylimits)
            options.ylimits(1) = min(min(y));
            options.ylimits(2) = max(max(y));
        end
        
        % Calculate offset
        y = y - options.i * (options.ylimits(2) * options.offset);
    end
end
end


% Check x-limits options
function options = plot_xlim(x, options)

% Determine y-limits
xmin = min(min(x));
xmax = max(max(x));

% Automatic x-limits
if isempty(options.xlimits)
    options.xlimits = [xmin, xmax];
    
% Manual x-limits
elseif strcmpi(options.xpermission, 'write')
    options.xlimits(1) = xmin;
    options.xlimits(2) = xmax;
end
end


% Check y-limits options
function options = plot_ylim(y, options)

% Determine y-limits
ymin = min(min(y));
ymax = max(max(y));

% Automatic y-limits
if isempty(options.ylimits)
    options.ylimits = [ymin, ymax];

% Manual y-limits
elseif strcmpi(options.ypermission, 'write')
    
    if ymin < options.ylimits(1)
        options.ylimits(1) = ymin;
    end
    if ymax > options.ylimits(2)
        options.ylimits(2) = ymax;
    end
end
end


% Update plot properties
function plot_update(options)

% Axes padding
padding.x = (options.xlimits(2) - options.xlimits(1)) * options.padding;
padding.y = (options.ylimits(2) - options.ylimits(1)) * options.padding;

% Axes limits
set(options.axes, 'xlim', [options.xlimits(1)-padding.x, options.xlimits(2)+padding.x]);
set(options.axes, 'ylim', [options.ylimits(1)-padding.y, options.ylimits(2)+padding.y]);

% Show legend
if strcmpi(options.legend, 'on')
    legend(options.axes, 'show', 'string', options.name);

% Hide legend
elseif strcmpi(options.legend, 'off')
    legend(options.axes, 'hide');
end
end


% Initialize axes
function options = plot_axes(obj, options, data)

% Axes properties
options.font.name = 'Avenir';
options.font.size = 14;
options.line.color = [0.23,0.23,0.23];
options.line.width = 1.25;
options.ticks.size = [0.007, 0.0075];

% Initialize figure
options.figure = figure(...
    'color', 'white',...
    'numbertitle', 'off',...
    'name', 'Chromatography Toolbox',...
    'units', 'normalized',....
    'position', obj.options.visualization.position,...
    'visible', 'off',...
    'paperpositionmode', 'auto');

% Determine color order
if ischar(options.colormap)
    try
        colors = colormap(options.colormap);
    catch
        options.colormap = 'parula';
        colors = colormap('parula');
    end
else
    colors = options.colormap;
end

% Check color group
if ~any(strcmpi(options.ions, {'tic', 'all'}))
    
    % Extracted ion chromatograms (selected)
    n = length(options.ions);
    
elseif strcmpi(options.ions, 'all')
    
    % Extracted ion chromatograms (all)
    n = length(data(options.samples(1)).mz);
else
    
    % Total ion chromatograms
    n = length(options.samples);
end

% Check amount of colors
if n < length(colors(:,1))
    
    % Use limited amount of colors
    c = round(length(colors)/n);
    options.colors = colors(1:c:end,:);
    
    if length(options.colors(:,1)) ~= n
        options.colors(1:length(options.colors) - n,:) = [];
    end
    
else
    % Use full amount of colors
    options.colors = colors;
end

% Initialize main axes
options.axes = axes(...
    'parent', options.figure,...
    'looseinset', [0.08, 0.1, 0.05, 0.05],...
    'fontsize', options.font.size-1,...
    'fontname', options.font.name,...
    'xcolor', options.line.color,...
    'ycolor', options.line.color,...
    'box', 'off',...
    'colororder', options.colors,...
    'color', 'none',...
    'linewidth', options.line.width,...
    'tickdir', 'out',...
    'ticklength', options.ticks.size,...
    'nextplot', 'add');

% Set x-axis label
options.xlabel = xlabel(...
    options.axes,...
    'Time (min)',...
    'fontsize', options.font.size,...
    'fontname', options.font.name);

% Overlaid y-axis label
if strcmpi(options.scale, 'normalized')
    options.ylabel = 'Intensity (%)';
else
    options.ylabel = 'Intensity';
end

% Stacked y-axis label
if strcmpi(options.layout, 'stacked') || options.offset ~= 0
    options.ylabel = [];
    set(options.axes, 'ytick', [], 'looseinset', [0.05, 0.1, 0.05, 0.05]);
end

% Set y-axis label
options.ylabel = ylabel(...
    options.axes,...
    options.ylabel,...
    'fontsize', options.font.size,...
    'fontname', options.font.name);

% Initialize empty axes
options.empty = axes(...
    'parent', options.figure,...
    'box','on',...
    'linewidth', options.line.width,...
    'color', 'none',...
    'xcolor', options.line.color,...
    'ycolor', options.line.color,...
    'xtick', [],...
    'ytick', [],...
    'position', get(options.axes, 'position'),...
    'nextplot', 'add');

box(options.empty, 'on');
box(options.axes, 'off');

% Link axes to allow zooming
linkaxes([options.axes, options.empty]);

% Version specific options
switch version('-release')
    
    case '2014b'
        
        % Resize callback
        set(options.figure, 'sizechangedfcn', @(varargin) set(options.empty, 'position', get(options.axes, 'position')));
        
        % Zoom callback
        set(zoom(options.figure), 'actionpostcallback', @(varargin) set(options.empty, 'position', get(options.axes, 'position')));

        % Axes overlap
        set(get(get(options.axes, 'yruler'),'axle'), 'visible', 'off');
        set(get(get(options.axes, 'xruler'),'axle'), 'visible', 'off');
        set(get(get(options.axes, 'ybaseline'),'axle'), 'visible', 'off');
        
    otherwise
        try
            % Resize callback
            set(options.figure, 'resizefcn', @(varargin) set(options.empty, 'position', get(options.axes, 'position')));
            
            % Zoom callback
            set(zoom(options.figure), 'actionpostcallback', @(varargin) set(options.empty, 'position', get(options.axes, 'position')));
        catch
        end
end
end


% Parse user input
function varargout = parse(varargin)

varargin = varargin{1};
nargin = length(varargin);

% Check input
if nargin < 1
    error('Not enough input arguments');
elseif isstruct(varargin{1})
    data = DataStructure('validate', varargin{1});
else
    error('Undefined input arguments of type ''data''');
end

% Check user input
input = @(x) find(strcmpi(varargin, x),1);

% Sample options
if ~isempty(input('samples'))
    samples = varargin{input('samples')+1};

    % Set keywords
    samples_all = {'default', 'all'};
        
    % Check for valid input
    if any(strcmpi(samples, samples_all))
        samples = 1:length(data);

    % Check input type
    elseif ~isnumeric(samples)
        
        % Check string input
        samples = str2double(samples);
        
        % Check for numeric input
        if ~any(isnan(samples))
            samples = round(samples);
        else
            samples = 1:length(data);
        end
    end
    
    % Check maximum input value
    if max(samples) > length(data)
        samples = samples(samples <= length(data));
    end
    
    % Check minimum input value
    if min(samples < 1)
        samples = samples(samples >= 1);
    end
    
    options.samples = samples;
    
else
    options.samples = 1:length(data);
end


% Ion options
if ~isempty(input('ions'))
    ions = varargin{input('ions')+1};
    
    % Set keywords
    ions_tic = {'default', 'tic', 'tics', 'total_ion_chromatograms'};
    ions_all = {'all', 'xic', 'xics', 'eic', 'eics', 'extracted_ion_chromatograms'};
    
    % Check for valid input
    if any(strcmpi(ions, ions_tic))
        options.ions = 'tic';
    
    elseif any(strcmpi(ions, ions_all))
        options.ions = 'all';

    elseif ~isnumeric(ions) && ~ischar(ions)
        options.ions = 'tic';
    else
        options.ions = ions;
    end

    % Check input range
    if isnumeric(options.ions)
        
        % Check maximum input value
        if any(max(options.ions) > cellfun(@length, {data(options.samples).mz}))
            options.ions = options.ions(options.ions <= min(cellfun(@length, {data(options.samples).mz})));
        end
        
        % Check minimum input value
        if min(options.ions) < 1
            options.ions = options.ions(options.ions >= 1);
        end
    end
    
else
    options.ions = 'tic';
end


% Baseline options
if ~isempty(input('baseline'))
    baseline = varargin{input('baseline')+1};
    
    % Set keywords
    baseline_on = {'on', 'show', 'display'};
    baseline_off = {'default', 'off', 'hide'};
    baseline_corrected = {'corrected', 'correct', 'subtract', 'subtracted'};
    
    % Check for valid input
    if any(strcmpi(baseline, baseline_off))
        options.baseline = 'off';
    
    % Check input type
    elseif any(strcmpi(baseline, baseline_on))
        options.baseline = 'on';
    
    elseif any(strcmpi(baseline, baseline_corrected))
        options.baseline = 'corrected';
        
    else
        options.baseline = 'off';
    end
    
else
    options.baseline = 'off';
end


% Peak options
if ~isempty(input('peaks'))
    peaks = varargin{input('peaks')+1};
    
    % Set keywords
    peaks_on = {'on', 'show', 'display'};
    peaks_off = {'default', 'off', 'hide'};
    peaks_residuals = {'residuals', 'residual', 'error', 'errors'};
    
    % Check for valid input
    if any(strcmpi(peaks, peaks_off))
        options.peaks = 'off';
    
    % Check input type
    elseif any(strcmpi(peaks, peaks_on))
        options.peaks = 'on';
    
    elseif any(strcmpi(peaks, peaks_residuals))
        options.peaks = 'corrected';
        
    else
        options.peaks = 'off';
    end
    
else
    options.peaks = 'off';
end


% Layout options
if ~isempty(input('layout'))
    layout = varargin{input('layout')+1};
    
    % Set keywords
    layout_stacked = {'stacked', 'stack', 'separate', 'separated'};
    layout_overlaid = {'default', 'overlaid', 'overlay', 'overlap', 'full'};
    
    % Check for valid input
    if any(strcmpi(layout, layout_stacked))
        options.layout = 'stacked';

    elseif any(strcmpi(layout, layout_overlaid))
        options.layout = 'overlaid';
        
    else
        options.layout = 'overlaid';
    end
    
else
    options.layout = 'overlaid';
end


% Scale options
if ~isempty(input('scale'))
    scale = varargin{input('scale')+1};
    
    % Set keywords
    scale_normalize = {'normalize', 'normalized', 'relative', 'separate', 'separated'};
    scale_full = {'default', 'full', 'all'};
    
    % Check for valid input
    if any(strcmpi(scale, scale_normalize))
        options.scale = 'normalized';
        
    elseif any(strcmpi(scale, scale_full))
        options.scale = 'full';
        
    else
        options.scale = 'full';
    end
    
else
    options.scale = 'full';
end


% Scaling options
if ~isempty(input('scaling'))
    scaling = varargin{input('scaling')+1};
    
    % Set keywords
    scaling_local = {'default', 'local', 'separated', 'separate', 'relative'};
    scaling_global = {'global', 'all', 'full'};
    
    % Check for valid input
    if any(strcmpi(scaling, scaling_local))
        options.scaling = 'local';
        
    elseif any(strcmpi(scaling, scaling_global))
        options.scaling = 'global';
        
    else
        options.scaling = 'local';
    end
    
else
    options.scaling = 'local';
end


% X-limits options
if ~isempty(input('xlim'))
    xlimits = varargin{input('xlim')+1};
    
    % Check for valid input
    if ~isnumeric(xlimits) || any(strcmpi(xlimits, {'default', 'auto'}))
        
        % Automatic x-limits
        options.xlimits = [];
        options.xpermission = 'write';
        
    % Check input length
    elseif xlimits(2) < xlimits(1) || length(xlimits) ~= 2;
    
        % Automatic x-limits
        options.xlimits = [];
        options.xpermission = 'write';
        
    else
        % Manual x-limits
        options.xlimits = xlimits;
        options.xpermission = 'read';    
    end

else
    options.xlimits = [];
    options.xpermission = 'write';
end


% Y-limits options
if ~isempty(input('ylim'))
    ylimits = varargin{input('ylim')+1};
    
    % Check user input
    if ~isnumeric(ylimits) || any(strcmpi(ylimits, {'default', 'auto'}))
        
        % Automatic y-limits
        options.ylimits = [];
        options.ypermission = 'write';
        
    % Check input length
    elseif ylimits(2) < ylimits(1) || length(ylimits) ~= 2
        
        % Automatic y-limits
        options.ylimits = [];
        options.ypermission = 'write';
    
    else 
        % Manual y-limits
        options.ylimits = ylimits;
        options.ypermission = 'read';
    end
    
else
    options.ylimits = [];
    options.ypermission = 'write';
end


% Axes padding
if ~isempty(input('padding'))
    padding = varargin{input('padding')+1};
    
    % Check for valid input
    if any(strcmpi(padding, {'default', 'on'}))
        options.padding = 0.05;
        
    elseif any(strcmpi(padding, {'off', 'none'}))
        options.padding = 0;

    % Check input range
    elseif padding(1) < 0 || padding(1) > 0.99
        options.padding = 0.05;
    else
        options.padding = padding(1);
    end
    
else
    options.padding = 0.05;
end
       

% Axes offset
if ~isempty(input('offset'))
    offset = varargin{input('offset')+1};
    
    % Check for valid input
    if any(strcmpi(offset, {'on'}))
        options.offset = 0.05;
        
    elseif any(strcmpi(offset, {'default', 'off', 'none'}))
        options.offset = 0;

    else
        options.offset = offset(1);
    end
    
else
    options.offset = 0;
end


% Linewidth options
if ~isempty(input('linewidth'))
    linewidth = varargin{input('linewidth')+1};
    
    % Check for valid input
    if strcmpi(linewidth, 'default') || ischar(linewidth)
        options.linewidth = 1.5;
        
    elseif linewidth <= 0
        options.linewidth = 1.5;

    else
        options.linewidth = linewidth(1);
    end
    
else
    options.linewidth = 1.5;
end


% Legend options
if ~isempty(input('legend'))
    legend = varargin{input('legend')+1};
    
    % Check for valid input
    if any(strcmpi(legend, {'default', 'off', 'hide'}))
        options.legend = 'off';
        
    elseif any(strcmpi(legend, {'on', 'show', 'display'}))
        options.legend = 'on';
        
    else
        options.legend = 'off';
    end
    
else
    options.legend = 'off';
end


% Colormap options
if ~isempty(input('colormap'))
    options.colormap = varargin{input('colormap')+1};
    
    % Check MATLAB version
    switch version('-release')
        case '2014b'
            default = 'parula';
            colormaps = {'parula', 'jet', 'hsv', 'hot', 'cool', 'spring', 'summer',...
                'autumn', 'winter', 'gray', 'bone', 'copper', 'pink', 'lines'};
        otherwise
            default = 'jet';
            colormaps = {'jet', 'hsv', 'hot', 'cool', 'spring', 'summer',...
                'autumn', 'winter', 'gray', 'bone', 'copper', 'pink', 'lines'};
    end

    % Check for valid input
    if isnumeric(options.colormap)
        
        % Check for RGB array
        if length(options.colormap(1,:)) ~= 3 || any(any(options.colormap > 1)) || any(any(options.colormap < 0))
            options.colormap = default;
        end
        
    % Check for default
    elseif strcmpi(options.colormap, 'default')
        options.colormap = default;
    
    % Check for colormap name
    elseif ~any(strcmpi(options.colormap, colormaps));
        options.colormap = default;
    end
    
else
    % Check MATLAB version
    switch version('-release')
        case '2014b'
            options.colormap = 'parula';
        otherwise
            options.colormap = 'jet';
    end
end


% Export options
if ~isempty(input('export'))
    export = varargin{input('export')+1};
    
    % Check for valid input
    if strcmpi(export, 'on')
        
        % Set default options
        options.export = {'chromatography_export', '-dpng', '-r300'};
    
    elseif any(strcmpi(export, {'default', 'off'}))
        options.export = [];
        
    % Check input type 
    elseif ~iscell(export)
        options.export = [];
    end
    
else
    options.export = [];
end

% Variables
options.name = {};

% Return input
varargout{1} = data;
varargout{2} = options;
end
