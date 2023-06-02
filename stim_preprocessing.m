clear all; close all; clc

%subject data parameters
subjID = 'S01';
subjFolder = 'C:\Users\chenli\Desktop\Beta_exp1\S01\';
fileNames = {'S01_1'; 'S01_2'; 'S01_3'; 'S01_4'; 'S01_5'; 'S01_6'; 'S01_7'; 'S01_8'; 'S01_9'; 'S01_10'};
cap = table2array(readtable("C:\Users\chenli\Desktop\cap.txt"));
outputDirectory = 'C:\Users\chenli\Desktop\Beta_exp1\indAvg_stim\';

channelsToAnalyze = [1:128];
reRef_channels = channelsToAnalyze;
downsample_pts = 2; 
numChannels = 1:128;

badCh_list_01 = [10 70 71 109];
badCh_list_02 = [10 70 71];
badCh_list_03 = [10 70 71];
badCh_list_04 = [10 70 71];
badCh_list_05 = [10 70 71];
badCh_list_06 = [105 108];
badCh_list_07 = [2 47 72 73];
badCh_list_08 = [7 34];
badCh_list_09 = [5];
badCh_list_10 = [11 69 70 71];
badCh_list_11 = [54 56 60 61 69 72];
badCh_list_12 = [31 32 34 72 73 74];
badCh_list_13 = [10 70 71];
badCh_list_14 = [32 34 92 122];
badCh_list_15 = [20 21 22 60 104 116];
badCh_list_16 = [10 46 57 70 71 97 103]; 

%frequency analyses parameters
bandpass_freq = [0.05 100];

%code parameters
stim_codes = [101 106 105 102]; %att stim loc 1, att stim loc 2, disatt stim loc 1, disatt stim loc 2
%att_codes = [101 102; 105 106]; %att loc 1, att loc 2
att_codes = [101 102; 105 106]; %att loc 1, att loc 2
cue_on_code = 20; %onset of the cue

%epoch timing parameters
preStim_time = -.5; %in seconds
postStim_time = .3; %in seconds
preStim_baseline = -.05; 
postStim_baseline = .03;
preCue_time = -.25; %in seconds
postCue_time = .7; %in seconds

interpolateData = 1; %variable to indicate whether data should be interpolated
artRec_th = 100;
rejCh_numberTrials_Threshold = 5;
interp_ch_Threshold = 20; %percentage of all trials
numChanInterp = 4;
badCh_list = badCh_list_01; 

%-------------------------------------------------------
%-------DON'T CHANGE ANYTHING BELOW HERE ---------------
for ss = 1 : length(fileNames) %substituting 'ss' variable to read each files 
    
    filename_bdf = [subjFolder fileNames{ss} '.bdf'];
    
    % *************************
    % read in EEG data bdf file
    % *************************
    %[EEG_data,numChan,labels,txt,fs] = eeg_read_bdf( filename_bdf,'all','n' );

    cfg = [];
    cfg.dataset = filename_bdf; %function that reads bdf??

    data = ft_preprocessing(cfg);
    EEG_data = cell2mat(data.trial);
    fs = data.fsample;
    labels = data.label;
    labels(end) = [];
    numChan = size(labels,1);

    % select only the trigger codes, not the battery and CMS status
    event = ft_read_event(filename_bdf);
    sel = find(strcmp({event.type}, 'STATUS'));
    event = event(sel);
    trigCodes = []; trigTimes = []; %both trigCodes and trigTimes are empty 
    trigCodes = [event.value]; %tell us what each trigger is based on the codes
    
    %extracting the channels to analyze
    EEG_data = EEG_data(channelsToAnalyze,:);

    %downsampling the EEG data (and making corresponding changes to the
    %FS-dependent variables)
    EEG_data = downsample(EEG_data',downsample_pts)';
    fs = fs ./ downsample_pts;
    trigTimes = round([event.sample]./downsample_pts);
    timeStim_vec = preStim_time:1/fs:postStim_time - 1/fs;
    timeCue_vec = preCue_time:1/fs:postCue_time - 1/fs;

    if ss == 1 
        %------------------
        %employing 50 hz noth filter
        % build notch filter
        nf = designfilt('bandstopiir','FilterOrder',4, ...
            'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
            'DesignMethod','butter','SampleRate', fs);
        %-----------------
        
        %designing a bandpass filter
        % and filtering the data
        BP_filter = designfilt('bandpassiir','FilterOrder',4, ...
            'HalfPowerFrequency1',bandpass_freq(1),'HalfPowerFrequency2',bandpass_freq(2), ...
            'SampleRate',fs);

        numCodes = length(find(trigCodes == stim_codes(1) | trigCodes == stim_codes(2) | ...
            trigCodes == stim_codes(3) | trigCodes == stim_codes(4)));
        all_epochs = zeros(size(numChannels,2),size(timeStim_vec,2), numCodes);
        typeTrials = cell(numCodes,2);
        numTrials = nan(length(numCodes),1);
   
        count = 1;

        %ERP_average
        ERP_data = zeros(size(all_epochs,1), size(all_epochs,2), length(stim_codes));
        countERP = zeros(length(stim_codes),1);


        %determining the location of the stimuli (around the 16 ring
        %configuration)

        fileToAnalyze = [subjFolder subjID '_1'];
        locDegrees = [];
        locDegrees = load([fileToAnalyze '.mat']).s.stim.locs;
        locDegrees = rad2deg(locDegrees);

    else
    end

    %filtering the data
    EEG_data = filtfilt(nf, EEG_data')';
    EEG_data = filtfilt(BP_filter, EEG_data')';

    %detrend
    numPts = fs * 4;
    startTime = 1:numPts:length(EEG_data)-fs*4;
    endTime = numPts:numPts:length(EEG_data);
    if endTime(end) ~=length(EEG_data)
        endTime(end) = length(EEG_data);
    end 

    for i = 1:length(startTime)
        EEG_data(startTime(i):endTime(i)) = detrend(EEG_data(startTime(i):endTime(i)));
    end 

    %-------------------------------

    %------------------------------
    %--------------------------------
    %--------------------------------
    %determining which channels to interpolate
    if interpolateData == 1
        chToInterp = badCh_list; 
%         chToInterp = [];
%         [chToInterp] = channelInterpolation_check(EEG_data,stim_codes,trigTimes,preStim_time,...
%             postStim_time,fs,timeStim_vec, artRec_th, length(numChannels),interp_ch_Threshold);

        %interpolating the data 
        distance = [];
        EEG_data_v2 = EEG_data;    
        time = 1:size(EEG_data_v2,2);
        for i = 1:length(chToInterp) %going through the bad channels
            c = chToInterp(i);
            distance = zeros(size(cap,1),2);
            distance(:,1) = cap(:,1);
            x = []; y = []; z = [];
            x = cap(c,2);
            y = cap(c,3);
            z = cap(c,4);
            distance(:,2) = sqrt((cap(:,2) - x).^2 + (cap(:,3) - y).^2  + (cap(:,4) - z).^2);
            distance = sortrows(distance,2);
            interpChannels = zeros(numChanInterp,1);
            if i < length(chToInterp)
               distance(ismember(distance(:,1), chToInterp(i+1:end)),:) = [];
            end 
            interpChannels = distance([2:(numChanInterp+1)],2);
            interpChannels = (sum(interpChannels) - interpChannels) / sum(interpChannels);
            EEG_data_v2(c,:) = sum(interpChannels .* EEG_data_v2(distance([2:numChanInterp+1],1),:),1);
        end
    else 
    end 

    %

    countTrialBlock = 1;

    %--------------
    %epoch stimuli only

    for nn = 1 : length(trigCodes)
        if trigTimes(nn) > 1
            curEpoch = []; tempCode = [];
            tempCode = find(stim_codes == trigCodes(nn));
            if ~isempty(tempCode)
    
                %epoch selection
                %-----------------
                curEpoch = EEG_data_v2(:,(trigTimes(nn) + round(preStim_time * fs)) : (trigTimes(nn) + round(postStim_time * fs)));
                if size(curEpoch,2) > length(timeStim_vec)
                    curEpoch(:,length(timeStim_vec)+1:end) = [];
                else
                end
                %-----------------
    
                %--------------------
                %baseline correction
                baselineVals = [];
                baselineVals = mean(curEpoch(:,find(timeStim_vec >= preStim_baseline & timeStim_vec < postStim_baseline)),2,'omitnan');
                curEpoch = curEpoch - repmat(baselineVals, 1, size(curEpoch,2));
    
                %----------------             
                %re-referencing data
                allChannels = numChannels;
                allChannels(badCh_list) = [];
                curEpoch = curEpoch - repmat(mean(curEpoch(allChannels,:), 1, 'omitnan'), size(curEpoch,1),1);

                %-----------------
                %flagging an art rej trial
                trialRej = 1;

                %doing rejection on electrodes that are considered 'good'
                curEpoch_goodChannels = [];
                curEpoch_goodChannels = curEpoch;
                curEpoch_goodChannels(badCh_list,:) = nan;

                if any(abs(curEpoch_goodChannels(:)) > artRec_th)
                    trialRej = -1;
                    [badCh,badPts] = find(abs(curEpoch) > artRec_th);
                    typeTrials{count,2} = unique(badCh);
                else
                    trialRej = 1;
                end
                
%                 if any(abs(curEpoch(:)) > artRec_th)
%                     trialRej = -1;
%                     [badCh,badPts] = find(abs(curEpoch) > artRec_th);
%                     typeTrials{count,2} = unique(badCh);
%                 else
%                     trialRej = 1;
%                 end
    
                %-----------------
    
                numTrials(count,1) = countTrialBlock;
                numTrials(count,2) = count;
                %numTrials(count) = nn;
                all_epochs(:,:,count) = curEpoch; 
                typeTrials{count,1} =  stim_codes(tempCode) * trialRej;
                count = count + 1;
                countTrialBlock = countTrialBlock + 1;
            else
            end
        end
   
    end
end

%----------------------
% deleting rejected trials
%----------------------
% deletedTrials = [];
typeTrialsCopy = typeTrials;
% for i = 1:size(typeTrialsCopy,1)
%     if length(typeTrialsCopy{i,2}) > rejCh_numberTrials_Threshold
%         deletedTrials = [deletedTrials; i];
%     else 
%     end 
% end

%typeTrialsCopy(deletedTrials,:) = [];

%--------
%averaging data for each condition
%--------
for cc = 1 : length(typeTrialsCopy)

    if typeTrialsCopy{cc,1} > 0 %then include (less than zeros i rejected trial)
        erpStim = [];
        erpStim = find(typeTrialsCopy{cc,1} == stim_codes);

        ERP_data(:,:,erpStim) = ERP_data(:,:,erpStim) + all_epochs(:,:,cc);
        countERP(erpStim) = countERP(erpStim) + 1;
    end
end

singleTrial = all_epochs;

%----------------------------
%re-referencing the data and making ind. subj averages
ERP_data_reRef = ERP_data;
for cc = 1 : length(stim_codes)
    ERP_data_reRef(:,:,cc) = ERP_data(:,:,cc) - ...
        repmat(squeeze(mean(ERP_data(reRef_channels,:,cc),1,'omitnan')),size(ERP_data,1),1);

    % ind. subject avg
    ERP_data_reRef(:,:,cc) = ERP_data_reRef(:,:,cc) ./ countERP(cc);
    ERP_data(:,:,cc) = ERP_data(:,:,cc)./ countERP(cc);
end
subj_File_Data = [];
subj_File_Data.ID = subjID;
subj_File_Data.trialTypes = typeTrialsCopy;
subj_File_Data.stimCodes = stim_codes;
subj_File_Data.preStimTime = preStim_time;
subj_File_Data.postStimTime = postStim_time;
subj_File_Data.timeVector = timeStim_vec;
subj_File_Data.samplingRate = fs;
subj_File_Data.RefChannels = reRef_channels;
subj_File_Data.locdeg = locDegrees; 
subj_File_Data.trialnumber = numTrials;

save([outputDirectory subjID '_indERP_data.mat'], 'ERP_data', 'ERP_data_reRef', 'singleTrial', 'countERP', 'subj_File_Data', '-v7.3');

%%
% %plotting to confirm interpolation works 
% figure;
% plot(1:length(all_epochs_v2{10}(128,:,84)),all_epochs_v2{10}(128,:,84));
% hold on; 
% plot(all_epoch_int{10}(128,:,84));
% plot(all_epochs_v2{10}(127,:,84));
% plot(all_epochs_v2{10}(126,:,84));
% plot(all_epochs_v2{10}(125,:,84));
% plot(all_epochs_v2{10}(7,:,84));
% hold off
% legend('old', 'new', '11', '12','10','9'); 
% hold on
% 
% %plotting an average of all the same stim code trials per channel 
% stimType = 101; 
% 
% for i = 1:10
%     for c = 1:128
%         stimTrials = find(abs(typeTrials(:,1)) == stimType);
%         plot(timeStim_vec, mean(all_epoch_int{i}(c,time,stimTrials),3));
%         title(num2str(c)); 
%         hold on
%         pause
%         close all
%     end
% end

