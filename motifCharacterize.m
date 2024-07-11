function [rst] = motifCharacterize(chValues, params)
    %% input parameter
    N = 0; % motif size (usually sqrt(# of bursts))
    thre = 0.01; % sec, similarity threshold
    mergethre = 0.9; % motif merge threshold
    validation = true; % flag for validation
    merge = true; % flag for merge
    assignthre = 3;
    display = true;

    inputidx = 1;
    while true
        if inputidx > length(params)
            break
        end
        
        switch params{inputidx}
            case 'N'
                inputidx = inputidx + 1;
                N = params{inputidx};
            case 'threshold'
                inputidx = inputidx + 1;
                thre = params{inputidx};
            case 'validation'
                inputidx = inputidx + 1;
                validation = params{inputidx};
            case 'merge'
                inputidx = inputidx + 1;
                merge = params{inputidx};
            case 'assign'
                inputidx = inputidx + 1;
                assignthre = params{inputidx};
            case 'display'
                inputidx = inputidx + 1;
                display = params{inputidx};
            otherwise
                error('MOTIFCHARACTERIZE: Input error');
        end
        inputidx = inputidx + 1;
    end
    
    rst.summary = [];
    rst.save = [];
    rst.visual = [];
    
    %% Activation matrix generation
    if ~ chValues.burstDetected % if burst detection
        burstDetect(chValues, {});
    end
    
    numActCh = nnz(chValues.active);
    burstnum = max(chValues.groups);
    bursts = struct('spktime', cell(burstnum, 1), 'chnum', cell(burstnum, 1));
    
    if burstnum <= 1
        return
    elseif burstnum > 400 % limit the number of burst to 400
        bursts = bursts(1:400);
        burstnum = 400;
    end
    
    for ii=1:burstnum
        idxtemp = (chValues.groups == ii);
        bursts(ii).spktime = chValues.timestamps(idxtemp);
        bursts(ii).chnum = chValues.chNums(idxtemp);
    end
    
    [actMats, isactive] = burst2actmat(bursts); % generation of activation matrix
    
    if N == 0
        N = floor(sqrt(burstnum));
    end
    
    %% Similarity matrix
    simMatOri = similarmat(actMats, thre); % calculate similarity matrix
    
    dist = 1 - simMatOri;
    dist = squareform(dist);
    z = linkage(dist, 'average');
    figure
    [~, ~, perm] = dendrogram(z, 0); % find permutation
    close
    simMat = simMatOri(perm, perm);
    actMats = actMats(perm); % reorder actMats according to hierarchical clustering
    isactive = isactive(perm);
    
    %% Motif finding
    [motifSet, motifActMat, motifPrcsn, motifCh] = motifSearch(actMats, isactive, N, thre);
    
    %% Motif validation
    simuldata = load('.\motif_simulation_total.mat');
    if validation
        % Check validity of each motif
        [valid, strd] = motifValidation(actMats, motifSet, motifActMat, motifPrcsn, motifCh, thre, simuldata, numActCh);
        
        motifSet = motifSet(valid);
        motifCh = motifCh(valid);
        strd = strd(valid);
    end
    
    
    %% Motif merging
    if merge
        % Check merge
        [ismerge, combs] = motifMerge(motifSet, motifActMat, motifPrcsn, motifCh, mergethre, numActCh);
        mergeto = 1:length(motifSet);
        
        for ii=1:length(ismerge)
            if ~ismerge(ii)
                continue
            end
            mergeto(combs(ii, 2)) = mergeto(combs(ii, 1));
        end
        
        motifSetRe = cell(max(mergeto), 1);
        for ii=1:max(mergeto)
            tempidx = find(mergeto == ii);
            for jj=1:length(tempidx)
                motifSetRe{ii} = cat(2, motifSetRe{ii}, motifSet{tempidx(jj)});
            end
        end
    else
        motifSetRe = motifSet;
    end
    
    % Recalculate activation matrix of motifs
    motifActMatRe = cell(length(motifSetRe), 1);
    motifPrcsnRe = cell(length(motifSetRe), 1);
    motifChRe = cell(length(motifSetRe), 1);
    for ii=1:length(motifSetRe)
        idxtemp = motifSetRe{ii};
        mattemp = zeros(size(actMats{1}));
        prcsntemp = zeros(size(actMats{1}));
        chtemp = true(length(isactive{1}), 1);
        
        for jj=1:length(idxtemp)
            mattemp = mattemp + actMats{idxtemp(jj)};
            prcsntemp = prcsntemp + actMats{idxtemp(jj)} .^ 2;
            chtemp = chtemp & isactive{idxtemp(jj)};
        end
        
        motifActMatRe{ii} = mattemp / length(idxtemp); % mean motif
        motifPrcsnRe{ii} = sqrt(prcsntemp / length(idxtemp) - motifActMatRe{ii} .^ 2); % std
        motifChRe{ii} = chtemp;
    end
    
    %% Burst allocation
    % Calcultate mean & std similarity of each motif
    meanSimilarity = zeros(length(motifSetRe), 1);
    stdSimilarity = zeros(length(motifSetRe), 1);
    for ii=1:length(motifSetRe)
        idxtemp = motifSetRe{ii};
        meansim = 0;
        stdsim = 0;
        for jj=1:length(idxtemp)
            simtemp = burstsim(motifActMatRe{ii}, actMats{idxtemp(jj)}, thre, 'chSelect', motifChRe{ii});
            meansim = meansim + simtemp;
            stdsim = stdsim + simtemp .^ 2;
        end
        meanSimilarity(ii) = meansim / length(idxtemp);
        stdSimilarity(ii) = sqrt(stdsim / length(idxtemp) - meanSimilarity(ii) .^ 2);
    end
    
    % find proper cluster 
    clusternum = zeros(burstnum, 1);
    for ii=1:burstnum
        simtemp = zeros(length(motifSetRe), 1);
        for jj=1:length(motifSetRe)
            if ismember(ii, motifSetRe{jj})
                clusternum(ii) = jj;
                break
            else
                simtemp(jj) = burstsim(motifActMatRe{jj}, actMats{ii}, thre, 'chSelect', motifChRe{jj});
            end
        end
        
        if clusternum(ii) > 0
            continue
        end
        
        [sim, I] = max(simtemp);
        if sim > meanSimilarity(I) - assignthre * stdSimilarity(I)
            clusternum(ii) = I;
        end
    end
    
    %% Make save data
    size(motifSetRe)
    size(motifActMatRe)
    motifstruct = struct('motifSet', motifSetRe, 'motifActMat', motifActMatRe, 'motifPrcsn', motifPrcsnRe, 'motifCh', motifChRe);
    motifsave = {simMat, simMatOri, burstnum, clusternum, motifstruct};

    %% Finish
    rst.summary = [];
    rst.save = motifsave;
    filenames = [];
    if display
        % similarty and dendrogram
        fig = figure('visible', 'off');
        imagesc(simMat)
        sgtitle(join(['Similarity (> ' num2str(thre) 'Hz)'],''));
    
        % similarity map
        subplot(3,2,1:4);
        imagesc(simMat)
        colormap(jet)
        colorbar
        pbaspect([1 1 1])
        xticklabels([])
        yticklabels([])
    
        % dendrogram
        simMatOri = simMat;
        subplot(3,2,5:6);
        dist = 1 - simMatOri;
        dist = squareform(dist);
        z = linkage(dist, 'average');
        h = dendrogram(z, 0);
        for ii=1:length(h)
            h(ii).Color = 'k';
        end
        xticklabels([])
        yticks([0 1])
        ylim([-0.01 1])
        
        nowstr = datestr(now, 'yymmdd-HHMMSS.FFF');
        filename = [nowstr '_similarity.png'];
        filename = fullfile(chValues.savepath, filename);
        filenames = [filenames filename];
        saveas(fig, filename);
        close(fig);

        % motif
        chWhole = [12, 13, 14, 16, 17, 21, 22, 23, 24, 25, 26, 27, 28, 31, 32, 33, 34, 35, 36, 37, 38, 41, 42, 43, 44, 45, 46, 47, 48, 51, 52, 53, 54, 55, 56, 57, 58, 61, 62, 63, 64, 65, 66, 67, 68, 71, 72, 73, 74, 75, 76, 77, 78, 82, 83, 84, 85, 86, 87];
        c = jet;
        c = flipud(c);
        maxtime = 0.03;

        for ii=1:length(motifstruct)
            actMat = motifstruct(ii).motifActMat;
            prcsn = motifstruct(ii).motifPrcsn;
            fstspk = firstSpikeRecon(actMat, prcsn);
            ts = fstspk(:, 1) * 1e3;
            prcsn = fstspk(:, 2) * 1e3;
            ts(ts > maxtime * 1e3) = NaN;
            ts(prcsn > maxtime / 2 * 1e3) = NaN;
            prcsn(ts > maxtime * 1e3) = NaN;
            prcsn(prcsn > maxtime / 2 * 1e3) = NaN;

            % subplot1 timestamp
            map = zeros(8) + maxtime + 0.1;
            map = map * 1e3;
        
            xs = zeros(length(chWhole), 1);
            ys = zeros(length(chWhole), 1);
            for jj=1:length(chWhole)
                x1 = floor(chWhole(jj) / 10);
                y1 = rem(chWhole(jj), 10);
                if ~isnan(ts(jj))
                    map(y1, x1) = ts(jj);
                end
                xs(jj) = x1;
                ys(jj) = y1;
            end

            fig = figure('visible', 'off');
            sgtitle(join(['Motif' num2str(ii)],''))
            subplot(311)
            imagesc(map)
            daspect([1 1 1])
            colormap(c);
            colorbar
            caxis([-0.005 maxtime] * 1e3)
            xticks([])
            yticks([])
            hold on
            [~, minidx] = min(ts);
            [~, maxidx] = max(ts);
            plot(xs(minidx), ys(minidx), 'xw', 'markersize', 10, 'linewidth', 3)
            plot(xs(maxidx), ys(maxidx), 'ow', 'markersize', 10, 'linewidth', 3)
            ylabel("spk time")

            % subplot2 order
            [ts, tsidx] = sort(ts);
            prcsn = prcsn(tsidx);
        
            map = zeros(8) + 60;
            map = map * 1e3;
        
            xs = zeros(length(chWhole), 1);
            ys = zeros(length(chWhole), 1);
            for jj=1:length(tsidx)
                x1 = floor(chWhole(tsidx(jj)) / 10);
                y1 = rem(chWhole(tsidx(jj)), 10);
                if ~isnan(ts(jj))
                    map(y1, x1) = jj;
                end
                xs(jj) = x1;
                ys(jj) = y1;
            end
            
            subplot(312)
            imagesc(map)
            daspect([1 1 1])
            colormap(c);
            colorbar
            caxis([0 30])
            xticks([])
            yticks([])
            hold on
            [~, minidx] = min(ts);
            [~, maxidx] = max(ts);
            plot(xs(minidx), ys(minidx), 'xw', 'markersize', 10, 'linewidth', 3)
            plot(xs(maxidx), ys(maxidx), 'ow', 'markersize', 10, 'linewidth', 3)
            ylabel("spk order")

            subplot(313)
            order = 1:59;
            errorbar(ts, order, prcsn, '.k', 'horizontal', 'linewidth', 2, 'markersize', 15)
            xlim([-.005 maxtime] * 1e3)
            ylim([0 30])
            ylabel("order"); xlabel("time [ms]")
            grid
            pbaspect([1 1 1])

            nowstr = datestr(now, 'yymmdd-HHMMSS.FFF');
            filename = [nowstr '_motif' num2str(ii) '.png'];
            filename = fullfile(chValues.savepath, filename);
            filenames = [filenames filename];
            saveas(fig, filename);
            close(fig);
        end
        
        rst.visual = filename;
    else
        rst.visual = '';
    end
end

