function [rst] = activeChs(chValues, params)
    %% inputsd
    inputidx = 1;
    threshold = 0.05; % Hz
    display = true;
    while true
        if inputidx > length(params)
            break
        end
        
        switch params{inputidx}
            case 'threshold'
                inputidx = inputidx + 1;
                threshold = params{inputidx};
            case 'display'
                inputidx = inputidx + 1;
                display = params{inputidx};
            otherwise
                error('ACTIVECHS: Input error');
        end
        inputidx = inputidx + 1;
    end
    
    %% calculation
    chs = chValues.getChs();
    active = false(length(chs), 1);
    
    for ii=1:length(chs)
        MFR = nnz(chValues.chNums == chs(ii)) / chValues.timespan;
        if MFR >= threshold
            active(ii) = true;
        end
    end
    chValues.active = active;
    chValues.activeChanneled = true;
    
    rst.summary = nnz(active);
    rst.save = chs(active);
    
    if display
        active8by8 = false(8, 8);
        for ii=1:length(chs)
            x = floor(chs(ii) / 10);
            y = rem(chs(ii), 10);
            if active(ii)
                active8by8(y, x) = true;
            end
        end

        fig = figure('visible', 'off');
        imagesc(active8by8);
        colormap jet
        colorbar
        daspect([1 1 1])
        xticklabels([])
        yticklabels([])
        title(['Active channels (> ' num2str(threshold) 'Hz)'])

        nowstr = datestr(now, 'yymmdd-HHMMSS.FFF');
        filename = [nowstr '_activechs.png'];
        filename = fullfile(chValues.savepath, filename);
        saveas(fig, filename);
        close(fig);
        rst.visual = filename;
    else
        rst.visual = '';
    end
end

