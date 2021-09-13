function FigsToPNG(directory,varargin)

argsin = varargin;

if isfolder(directory)
    cd(directory)
    files = dir('*.fig');
    for i = 1:length(files)
        fig = openfig(files(i).name);
        name = erase(files(i).name,'.fig');
        set(gcf, 'Color', 'w');
        export_fig([name '.png'],argsin{:})
        %saveas(fig,[name '.jpg']);
        close(fig)
    end
elseif isfile(directory)
    fig = openfig(directory);
    name = erase(directory,'.fig');
    set(gcf, 'Color', 'w');
    export_fig([name '.png'],argsin{:})
    %saveas(fig,[name '.jpg']);
end