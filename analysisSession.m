classdef analysisSession < handle
    properties
        filelist
        filenum
        
        config
        analysisConfig
        analysisnum
        analysisname
        analysisinputs
        
        analysisSummary
        analysisSave
        analysisVisual
        txtPrint
        
        intervalConfig
        
        savepath
    end
    
    methods
        function obj = analysisSession(filelist)
            obj.filelist = filelist;
            obj.filenum = length(filelist);
            
            obj.config = readtable('analysis_base.xlsx');
            obj.intervalConfig = {};
            analysisSession.loadLibrary();
        end
        
        function loadAnalysisConfig(obj, configFile)
            if isempty(configFile)
                configFile = 'analysis_config.xlsx';
            end
            if ~ischar(configFile) || exist(configFile, 'file') ~= 2
                exist(configFile, 'file')
                error('Analysis config error');
            end
            
            opts = detectImportOptions(configFile);
            opts.DataRange = 'A1';
            obj.analysisConfig = readcell(configFile, opts);
            
            obj.analysisnum = size(obj.analysisConfig, 1);
            obj.analysisname = obj.analysisConfig(:, 1);
            
            obj.analysisinputs = cell(obj.analysisnum, 1); % analysis inputs
            for ii=1:obj.analysisnum
                input = cell(size(obj.analysisConfig, 2) - 1, 1);
                for jj=2:size(obj.analysisConfig, 2)
                    if ismissing(obj.analysisConfig{ii, jj})
                        input(jj - 1:end) = [];
                        obj.analysisinputs{ii} = input;
                        break
                    else
                        input{jj - 1} = obj.analysisConfig{ii, jj};
                    end
                    obj.analysisinputs{ii} = input;
                end
            end
            
            obj.analysisSummary = false(obj.analysisnum, 1);
            obj.analysisSave = false(obj.analysisnum, 1);
            obj.analysisVisual = false(obj.analysisnum, 1);
            obj.txtPrint = strings(obj.analysisnum, 1);
            for ii=1:obj.analysisnum
                methodidx = obj.findMethod(obj.analysisname{ii});
                obj.analysisSummary(ii) = obj.config.Summary(methodidx);
                obj.analysisSave(ii) = obj.config.Save(methodidx);
                obj.analysisVisual(ii) = obj.config.Visualize(methodidx);
                obj.txtPrint(ii) = obj.config.txtPrint(methodidx);
            end
        end
        
        function analyze(obj)
            if isempty(obj.filelist)
                error('ANALYZE: There are no files to analyze')
            end
            nowstr = datestr(now, 'yymmdd-HHMM');
            obj.savepath = fullfile('.', nowstr);
            [status, msg] = mkdir(obj.savepath);
            if status ~= 1
                error(['ANALYZE: folder making failed, ' msg])
            end
            total = cell(obj.filenum, obj.analysisnum);
            metas = cell(obj.filenum, 1);
            for ii=1:obj.filenum
                filepath = fullfile(obj.filelist(ii).folder, obj.filelist(ii).name);
                myvalue = obj.fileLoad(filepath);
                myvalue.savepath = obj.savepath;

                for jj=1:length(obj.intervalConfig)
                    intvConfig = obj.intervalConfig{jj};
                    myvalue.buildInterval(intvConfig.name, intvConfig.manner, intvConfig.input1,  intvConfig.input2);
                end
                metas{ii} = myvalue.meta;
                for jj=1:obj.analysisnum
                    analysisName = obj.analysisname{jj};
                    inputs = obj.analysisinputs{jj};
                    total{ii, jj} = obj.connectMethod(analysisName, myvalue, inputs);
                end
            end
            
            % Reorganize analysis names
            anname = cell(obj.analysisnum, 1);
            anname{1} = obj.analysisname{1};
            for ii=2:length(obj.analysisname)
                dupl = sum(strcmp(obj.analysisname(1:ii - 1), obj.analysisname{ii}));
                if dupl > 0
                    anname{ii} = [obj.analysisname{ii} num2str(dupl)];
                else
                    anname{ii} = obj.analysisname{ii};
                end
            end
            obj.analysisname = anname;
            
            % Need to implement data save, summary, visualization
            obj.savedata(total, metas);
            obj.summary(total, metas);
            obj.summary_(total);
            obj.visualization(total);
        end
        
        function myvalue = fileLoad(obj, filepath)
            % file load
            if endsWith(filepath, 'mcd') % spike data
                myvalue = chValues(filepath);
            elseif endsWith(filepath, 'h5')
                myvalue = chValues(filepath);
            elseif endsWith(filepath, 'gdf')
                myvalue = chValues(filepath);
            elseif endsWith(filepath, 'corr') % correlation data
                myvalue = corrValues(filepath);
            end
        end
        
        function rst = connectMethod(obj, analysisName, myvalue, inputs)
            % Practical execution of methods
            methodName = obj.findMethodName(analysisName);
            
%             inputStr = ''; % Concat input args as input string
%             for ii=1:length(inputs)
%                 if ischar(inputs{ii})
%                     inputStr = strcat(inputStr, inputs{ii});
%                 elseif isnumeric(inputs{ii})
%                     inputStr = strcat(inputStr, num2str(inputs{ii}));
%                 end
%                 
%                 if ii ~= length(inputs)
%                     inputStr = strcat(inputStr, ', ');
%                 end
%             end
            
            rst = eval([methodName, '(myvalue, inputs);']);
        end
        
        function methodidx = findMethod(obj, analysisName)
            names = obj.config.AnalysisName;
            methodidx = strcmp(names, analysisName);
            
            if sum(methodidx) == 0
                error(['FINDMETHOD: check analysis name: ' analysisName])
            end
        end
        
        function methodName = findMethodName(obj, analysisName)
            methodidx = obj.findMethod(analysisName);
            methods = obj.config.MethodName;
            methodName = methods{methodidx};
        end
        
        function savedata(obj, total, metas)
            % implemetation
            datatosave = cell(size(total));
            for ii=1:obj.filenum
                for jj=1:obj.analysisnum
                    datatosave{ii, jj} = total{ii, jj}.save;
                end
            end
            
            files = obj.filelist;
            analysis = obj.analysisname;
            inputs = obj.analysisinputs;
            
            filename = analysisSession.savenameGen(obj.savepath);
            save(filename, 'datatosave', 'files', 'analysis', 'inputs', 'metas');
        end

        function summary_(obj, total)
            for jj=1:obj.filenum
                filenameOri = obj.filelist(jj).name;
                filename = analysisSession.summarytxtnameGen(obj.savepath,filenameOri);
    
                summaryStr = filenameOri;
                summaryStr = strcat(summaryStr,"\n\n");
                
                for ii=1:obj.analysisnum
                    tempcell = [];
                    if ~obj.analysisSummary(ii)
                        continue
                    end
                    sumsum = total{jj, ii}.summary;
                    if isempty(sumsum)
                        continue
                    end

                    summaryStr = strcat(summaryStr,obj.txtPrint(ii));
                    summaryStr = strcat(summaryStr,"\n");
                    summaryStr = strcat(summaryStr,num2str(sumsum));
                    summaryStr = strcat(summaryStr,"\n\n");
                end
                
                % write
                fileID = fopen(filename,'w');
                fprintf(fileID,summaryStr);
                fclose(fileID);
            end
        end
        
        function summary(obj, total, metas)
            % implementation, addvars
            summaryTable = table();
            
            metanames = fieldnames(metas{1});
            tempcell = cell(obj.filenum, length(metanames));
            
            for ii=1:obj.filenum
                for jj=1:length(metanames)
                    tempcell{ii, jj} = metas{ii}.(metanames{jj});
                end
            end

            summaryTable = addvars(summaryTable, tempcell);
            summaryTable = splitvars(summaryTable);
            summarylabels = cell(obj.analysisnum, 1);
            tempidx = 1;
            for ii=1:obj.analysisnum
                tempcell = [];
                if ~obj.analysisSummary(ii)
                    continue
                end
                for jj=1:obj.filenum
                    sumsum = total{jj, ii}.summary;
                    if isempty(sumsum)
                        continue
                    end
                    tempcell = cat(1, tempcell, sumsum);
                end
                if isempty(tempcell)
                    continue
                end
                summaryTable = addvars(summaryTable, tempcell);
                summarylabels{tempidx} = obj.analysisname{ii};
                tempidx = tempidx + 1;
            end
            summarylabels(tempidx:end) = [];
            summaryTable.Properties.VariableNames = cat(1, metanames, summarylabels);
            
            save("summaryTable.mat","summaryTable");

            conditionTable = table(obj.analysisConfig);
            
            filename = analysisSession.summarynameGen(obj.savepath);
            writetable(summaryTable, filename, 'Sheet', 'Summary');
            writetable(conditionTable, filename, 'Sheet', 'Condition');
        end
        
        function visualization(obj, total)
            import mlreportgen.ppt.*
            
            filename = analysisSession.visualnameGen(obj.savepath);
            ppt = Presentation(filename);
            open(ppt);
            
            titleslide = add(ppt, 'Title Slide');
            tstr = [obj.filelist(1).folder '\nCharacterization'];
            t = Text(tstr);
            t.FontSize = '40pt';
            p = Paragraph;
            append(p, t);
            replace(titleslide, 'Title', p);
            nowstr = datestr(now, 'yymmdd-HHMM');
            replace(titleslide, 'Subtitle', nowstr);
            
            for ii=1:obj.filenum
                contslide = add(ppt, 'Title Only');
                t = Text(obj.filelist(ii).name);
                t.FontSize = '24pt';
                p = Paragraph;
                append(p, t);
                replace(contslide, 'Title', p);
                
                coordi = 0;
                for jj=1:obj.analysisnum
                    if ~obj.analysisVisual(jj)
                        continue
                    end
                    pic = total{ii, jj}.visual;
                    if isempty(pic)
                        continue
                    end
                    rstpic = Picture(pic);
                    rstpic.Height = '6cm';
                    rstpic.Width = '8cm';
                    
                    x = rem(coordi, 4);
                    y = floor(coordi / 4);
                    rstpic.X = [num2str(x * 8 + 0.5) 'cm'];
                    rstpic.Y = [num2str(y * 6 + 5) 'cm'];
                    
                    add(contslide, rstpic);
                    coordi = coordi + 1;
                end
            end
            
            close(ppt);
        end
        
        function buildInterval(obj, name, manner, input1, input2)
            intvConfig.name = name;
            intvConfig.manner = manner;
            intvConfig.input1 = input1;
            intvConfig.input2 = input2;
            obj.intervalConfig{end + 1} = intvConfig;
        end
    end
    
    methods (Static)
        function loadLibrary()
%             libraryPath = [".\mcd\nsMCDlibrary64.dll", ...
%             "..\mcd\nsMCDlibrary64.dll", ...
%             "C:\Program Files\MATLAB\R2020a\toolbox\mcd\nsMCDlibrary64.dll", ...
%             "D:\Program Files\MATLAB\R2020a\toolbox\mcd\nsMCDlibrary64.dll", ...
%             "C:\Program Files\MATLAB\R2021a\toolbox\mcd\nsMCDlibrary64.dll", ...
%             "D:\Program Files\MATLAB\R2021a\toolbox\mcd\nsMCDlibrary64.dll"];
%     
%             stat = 0;
%             for ii=1:length(libraryPath)
%                 if exist(libraryPath(ii), 'file') == 2
%                     stat = 1;
%                     libpath = libraryPath(ii);
%                     break
%                 end
%             end
            libfile = 'nsMCDlibrary64.dll';
            libpath = which(libfile);

            if isempty(libpath)
                error('LOADLIBRARY: library not found')
            end

            ns = ns_SetLibrary(char(libpath));
            if ns ~= 0
                error('LOADLIBRARY: library load failed');
            end
        end
        
        function filename = savenameGen(savepath)
            nowstr = datestr(now, 'yymmdd-HHMM');
            filename = fullfile(savepath, [nowstr '.mat']);
        end
        
        function filename = summarynameGen(savepath)
            nowstr = datestr(now, 'yymmdd-HHMM');
            filename = fullfile(savepath, [nowstr '.xlsx']);
        end

        function filename = summarytxtnameGen(savepath,filenameOri)
            nowstr = datestr(now, 'yymmdd-HHMM');
            filename = fullfile(savepath, [nowstr '_' filenameOri '.txt']);
        end
        
        function filename = visualnameGen(savepath)
            nowstr = datestr(now, 'yymmdd-HHMM');
            filename = fullfile(savepath, [nowstr '.pptx']);
        end
    end
end

