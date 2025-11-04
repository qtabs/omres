function omres(blockNumber)

    % omRes session including functional localiser and oddball task
    %   blockNumber = 0 or 1 to start from the beginning (default)
    %               = 2, 3, or 4 to start from blocks 2, 3, or 4
    %               = -1 to resume from last saved cached data
    %               = 'soundcheck' (or any other string) for sound check

    if nargin < 1, blockNumber = 0; end;

    pars = getParameters();

    if ischar(blockNumber)
        if strcmp(blockNumber, 'findport')
            findPort(pars);
        else
            soundCheck(pars);
        end
        return
    end

    mess    = getMessages(pars);
    session = initialiseExp(blockNumber, pars, mess);
    
    fprintf('\nRunning subject "%s", ',  session.subject)
    fprintf('session #%d, ',           session.sessN)
    fprintf('starting at block #%d. ', session.currentN)
    input('Press ENTER to start...')
    fprintf('\n')

    handles = initialiseHandles(pars);
    %handles = initialiseHandlesLight(pars);

    try
    	session = runSession(session, pars, mess, handles);
        closeInterfaces();
        delete(pars.cache);
    catch e
    	closeInterfaces();
    	rethrow(e);
    end

end



function session = runSession(session, pars, mess, handles)

    n0 = session.currentN;

    if ismember('localiser', pars.task(pars.scheme))
        sounds = readSounds(pars);
    end

    for n = n0:length(session.runs)

        session.currentN = n;
        if session.currentN > 1
            s = session.score;
        end
        save(pars.cache, '-struct', 'session');

        fprintf('\nReady! Next in line: %s ', pars.task{pars.scheme(n)});
        fprintf('(Block #%d). ', n);
        fprintf('Press space to start...\n');
        waitForSpace(handles);

        [handles.logF, handles.messF] = openLogfiles(session, n, pars);
        handles.absOnset = waitForScanner(handles, mess);

        if strcmp(pars.task{pars.scheme(n)},'omissions')
            s(n) = runOmissions(session.runs{n},n,'tones',pars,mess,handles);
        elseif strcmp(pars.task{pars.scheme(n)},'localiser')
            s(n) = runLocaliser(session.runs{n},sounds,pars,handles);
        end

        session.score = s;

        try, fclose(handles.logF); fclose(handles.messF); end;

        if n == length(session.runs)
            displayText(handles, mess.endExp, false);
            perfStr = [sprintf('%.2f, ', s), repmat(char(8), [1, 2])];
            score   = mean(s);
            fprintf('\n[Perf = [%s]; Avg perf = %.2f]\n\n', perfStr, score);
            consciousWait(60, handles)
        else
            displayText(handles, mess.waitRun, false);
        end

    end

end



function score = runOmissions(blck, blckIx, stim, pars, mess, handles)

    performance = [];

    for i = 1:size(blck, 1)

        if blck(i, 2) == -1
            fprintf('[%s %2d/%d] (null)\n', mess.trial, i, size(blck, 1));
        elseif blck(i, 2) == 0
            fprintf('[%s %2d/%d] (no omission)\n', mess.trial, i, size(blck, 1));
        else
            fprintf('[%s %2d/%d]\n', mess.trial, i, size(blck, 1));
        end

        correct = playTrain(blck, i, pars, mess, handles, stim);
        
        if not(isempty(correct))
            performance = [performance, correct];
        end

    end

    score = mean(performance);

end



function score = runLocaliser(blck, sounds, pars, handles)

    performance = [];

    for j = 1:length(blck)

        fprintf('Block %d/%d\n',  j, length(blck));
        consciousWait(0.5, handles);
        displayText(handles, '+')

        for i = 1:length(blck{j})  
            event = blck{j}(:, i);   
            [onset, tKey, result, ID] = playSound(event, sounds, pars, handles);
            writeSoundToLog(handles, onset, 1.0, tKey, result, ID);

            if strcmp(result, 'detection')
                performance = [performance, 1];
            elseif strcmp(result, 'miss') || strcmp(result, 'falseAlarm')
                performance = [performance, 0];
            end
        end

        IBI = min(pars.loc.maxIBI, max(pars.loc.minIBI, pars.loc.avgIBI+randn()));
        consciousWait(1, handles);
        displayText(handles, '+')
        consciousWait(max(0, IBI - 1), handles)
    end

    score = mean(performance);

end



function [tSound, tKey, result, ID] = playSound(event, sounds, pars, handles);

    if event(3) == 0
        ID = 'silence';
        fprintf(' --- [%s]\n', ID) 
        tSound = GetSecs() - handles.absOnset;
        tKey = -1;
        result = '';
        consciousWait(1 + pars.loc.ISI, handles);
    else
        ID = sounds{event(1)}.label;
        s  = pars.natsoundGain * sounds{event(1)}.waveform;
        fprintf(' --- [%s]\n', ID)  
        PsychPortAudio('FillBuffer', handles.audio, s);       
        PsychPortAudio('Start', handles.audio); 
        
        [tKey, pressed] = funLocFeedback(pars, handles);
        if pressed
            tKey = tKey - handles.absOnset;
            if event(2)
                result = 'detection';
                displayText(handles, '+green')
            else
                result = 'falseAlarm';
                displayText(handles, '+red')
            end
        else
            tKey = -1;
            if event(2)
                result = 'miss';
                displayText(handles, '+red')
            else
                displayText(handles, '+')
                result = '';
            end
        end
        
        tSound = PsychPortAudio('Stop', handles.audio, 1, 1);
        tSound = tSound - handles.absOnset;
    end

end



function performance = updatePerformance(performance, result)

    %performance = [trueDetections, missedRep, falseAlarms]

    if nargin == 1
        performance = performance(1) / (sum(performance) + eps);
    elseif strcmp(result, 'detection')
        performance(1) = performance(1) + 1;
    elseif strcmp(result, 'miss')
        performance(2) = performance(2) + 1;
    elseif strcmp(result, 'falseAlarm')
        performance(3) = performance(3) + 1;
    end

end



function correct = playTrain(blck, i, pars, mess, handles, stim)

    %blck
    if blck(i, 2) == -1 % Null event
        duration = pars.ob.len * (sum(pars.ob.dur) + pars.ob.ISI);
        duration = duration + pars.ob.ISI;
        onsets   = GetSecs() - handles.absOnset;
        consciousWait(duration, handles);
        a = 0; rt = 0; correct = [];
        comb = sprintf('null');
    else
        sStn = createSinusoid(pars.ob.tones(blck(i, 1)), pars);
        s   = [];
        for j = 1:pars.ob.len
            if j == blck(i, 2)
                s = [s, 0 * sStn.s];
            else 
                s = [s, sStn.s];       
            end
            if j < pars.ob.len
                s = [s, zeros([1, pars.ob.ISI * pars.fs])];
            end
        end

        PsychPortAudio('FillBuffer', handles.audio, s);
        PsychPortAudio('Start', handles.audio);

        [a, correct, rt] = omissionsFeedback(blck(i,2), pars, mess, handles);
    
        startTime = PsychPortAudio('Stop', handles.audio, 1, 1);

        onsets = (0:1:(pars.ob.len - 1)) * (sum(pars.ob.dur) + pars.ob.ISI);
        onsets = onsets + startTime - handles.absOnset;

        comb = sprintf('std%dHz', round(pars.ob.tones(blck(i, 1))));
    end
    
    writeTrainToLog(handles, onsets, blck(i, :), comb, a, rt, pars);  
    
    if blck(i, 2) > -1
        consciousWait(computeITI(blck, i, pars, rt), handles);
    end

end



function [answer, correct, rt] = omissionsFeedback(loc, pars, mess, handles)

    if loc == 0 % If no omission, participant will know after 6 stds
        timeout = pars.ob.len*(pars.ob.ISI+sum(pars.ob.dur))-pars.ob.ISI;
        answerTimeWatch = GetSecs();
    else        
        devOnset = loc*(pars.ob.ISI+sum(pars.ob.dur))-pars.ob.ISI;
        timeout  = devOnset + pars.ob.maxRT;
        answerTimeWatch = GetSecs() + devOnset;
    end
   
    [keyPressed, keyPressTime] = waitForKey(handles, timeout);
    
    if not(isempty(keyPressed))
        keyTime = keyPressTime - handles.absOnset;
        %keyEvent = ['buttonPress-', keyPressed(find(~isspace(keyPressed)))];
        keyEvent = 'buttonPress';
        logInput = sprintf('%.4f\t%.4f\t%s\n', keyTime, 0.1, keyEvent);
        fprintf(handles.logF, logInput);
        rt = keyPressTime - answerTimeWatch;
    else
        rt = [];
    end

    if isempty(keyPressed)
        if loc == 0
            displayText(handles, '+green');
            answer = 0;
            correct = true;
        else
            displayText(handles, '+red');
            answer = 0;
            correct = false;
        end
        consciousWait(pars.ob.feedbackT, handles);
    elseif length(keyPressed) > 1
        displayText(handles, '+red');
        consciousWait(pars.ob.feedbackT, handles);
        answer = -1;
        correct = false;
    else
        answer = keyToButton(keyPressed, handles) + min(pars.ob.devLocs) - 1;
        if answer == loc
            displayText(handles, '+green');
            correct = true;
            consciousWait(pars.ob.feedbackT, handles);
        else
            displayText(handles, '+red');
            consciousWait(pars.ob.feedbackT, handles);
            correct = false;
        end
    end

    displayText(handles, '+');

end



function [pressTime, pressed] = funLocFeedback(pars, handles) 

    t0 = GetSecs();
    pressed   = false;
    pressTime = 0;
    waitTime  = pars.loc.maxRT;

    [~, keyPressTime, keyCode] = KbCheck();
    %keyPressTime = GetSecs();

    while (sum(keyCode(handles.keys.button))==0) && (GetSecs()-t0 < waitTime)
        [~, keyPressTime, keyCode] = KbCheck();
        
        if (keyCode(handles.keys.escape))
            terminateExecution();     
        elseif sum(keyCode(handles.keys.button)) > 0
            if not(pressed)
                pressTime = keyPressTime;
                pressed = true;
            end
        end
    end

    consciousWait(GetSecs()-t0 - waitTime, handles);

end



function stim = createSinusoid(f, pars)
   
    N = pars.ob.dur * pars.fs;
    g = 1/loudWeight(f, pars.phons);

    fTrans = linspace(f(1), f(end), N(2)); 
    stim.f = [f(1) * ones([1, N(1)]), fTrans, f(end) * ones([1, N(3)])];
    stim.s = g * pars.beepGain * sin(2 * pi * cumsum(stim.f) / pars.fs);
    stim = ramps(stim, pars.ob.dur(1), pars.fs);

end



function stim = ramps(stim, T, fs)

    if nargin < 3, fs = 48000; end;
    if nargin < 2, T = 0.05; end;

    N = T * fs;
    
    hamm = hamming(2 * N);
    modd = ones(size(stim.s));

    modd(1:N) = hamm(1:N);
    modd((length(stim.s) - N + 1):end) = hamm((N+1):end);
    stim.s = stim.s .* modd;

end



function runs = generateRuns(pars)

    if ismember('omissions', pars.task(pars.scheme))
        omissionRuns  = generateOmissionRuns(sum(pars.scheme==1), pars);
    end
	
	tix = 1; six = 1;
	for i = 1:length(pars.scheme)
        if strcmp(pars.task{pars.scheme(i)}, 'omissions')
            runs{i} = omissionRuns{tix};
            tix = tix + 1;
        elseif strcmp(pars.task{pars.scheme(i)}, 'localiser')
            runs{i} = generateLocaliserRun(pars);
        end
	end

end



function blcks = generateOmissionRuns(nOfRuns, pars)

    combs = repmat(pars.ob.standards, [length(pars.ob.devLocs), 1]);
    combs = combs(:);
    combs = [combs, repmat(pars.ob.devLocs', [length(pars.ob.standards), 1])];

    totalNoOmissions = round(nOfRuns * size(combs, 1) * pars.ob.probNoOm);
    noOmsPerStandard = totalNoOmissions / length(pars.ob.standards);
    nNoOmissions  = zeros([nOfRuns], 1);
    noCombs = repmat(pars.ob.standards', [noOmsPerStandard, 1]);
    noCombs = [noCombs, repmat(0, [length(noCombs), 1])];
    noCombs = noCombs(randperm(size(noCombs, 1)), :);

    noOmTrials = {[], [], [], []}; k = 0;
    while k < totalNoOmissions
        noOmTrials{mod(k,4)+1} = [noOmTrials{mod(k,4)+1}; noCombs(k+1, :)];
        k = k + 1;
    end

    for i = 1:nOfRuns
        nOfNullEvs = round((size(combs,1) + length(noOmTrials{i})) * pars.ob.null);
        nOfWeirdEvents = nOfNullEvs + length(noOmTrials{i});
        nOfEvents = size(combs, 1) + nOfWeirdEvents;
        
        weirdTrials = 2 * randperm(floor(nOfEvents / 2), nOfWeirdEvents);
        nullEvents = weirdTrials(1:nOfNullEvs);
        noOmEvents = weirdTrials((nOfNullEvs+1):end);
        
        blcks{i} = combs(randperm(size(combs, 1)), :);
        k = 1;
        for j = sort(weirdTrials)
            blcks{i}((j+1):(size(blcks{i}, 1)+1), :) = blcks{i}(j:end, :);
            if ismember(j, nullEvents)
                blcks{i}(j, :) = [-1, -1];
            elseif ismember(j, noOmEvents)
                blcks{i}(j, :) = noOmTrials{i}(k, :);
                k = k + 1;
            end
        end
    end

end



function blcks = generateLocaliserRun(pars)


    % 1. Gather all the sounds in the stimuli path
    soundDir  = dir([pars.loc.stimPath '/*s3_*.wav']);

    % 2. Create a vector containing a number for each sound sIx = [1 ... 84] 
    %    We will use these numbers to let "playStim" know which sound we want
    %    to play each time
    sIx = 1:length(soundDir);

    % 3. Create a pool of all the sounds (sound IDs) we will need in total 
    sIxPool = [];
    while length(sIxPool) < (pars.loc.Nblocks * pars.loc.blockLen)
        sIx = sIx(randperm(length(sIx))); % Shuffles sIX to randomise the order
        sIxPool = [sIxPool, sIx];         % Adds the IDs to the pool
    end

    % 4. Split the pool in pars.Nblocs blocks with pars.blockLen elements
    %    Note that sIxPool might be longer than the total amount of sounds
    %    we need, that's why we shuffle sIx all the time, to randomise which
    %    sounds appear less times in each run
    for i = 1:pars.loc.Nblocks
        ix0 = 1 + (i - 1) * pars.loc.blockLen;
        ix1 = i * pars.loc.blockLen;
        soundBlcks{i}(1, :) = sIxPool(ix0:ix1);
        soundBlcks{i}(3, :) = 1; % No silence
    end

    % Now we have an array of blocks, "soundBlcks". The n^th element,
    % soundBlcks{n}, contains the IDs of the sounds we will play in that block

    % 5. Insert repetitions with a probability repProb
    for i = 1:pars.loc.Nblocks
        for j = 2:pars.loc.blockLen
            if rand() < pars.loc.repProb
                soundBlcks{i}(1, j) = soundBlcks{i}(1, j - 1);
                soundBlcks{i}(2, j) = 1;
            else % if there was a repetition already we need to tackle it
                if soundBlcks{i}(1, j) == soundBlcks{i}(1, j - 1)
                   soundBlcks{i}(1,j) = mod(soundBlcks{i}(1,j),length(sIx))+1;
                end
                soundBlcks{i}(2, j) = 0;
            end
        end
    end

    % 6. We check that there are indeed as many repetitions as we wanted:
    %for i = 1:pars.loc.Nblocks
    %    rep(i) = sum((soundBlcks{i}(2:end) - ...
    %                    soundBlcks{i}(1:(pars.loc.blockLen-1))) == 0);
    %end
    %sum(rep) / (pars.loc.Nblocks * (pars.loc.blockLen - 1))

    % 7. Generate blocks by alternating soundblocks with silent blocks
    silentBlck = repmat([0; 0; 0], [1, pars.loc.blockLen]);
    for i = 1:(2 * pars.loc.Nblocks)
        if mod(i, 2) == 1
            blcks{i} = soundBlcks{(i+1) / 2};
        else
            blcks{i} = silentBlck;
        end
    end

end



function session = initialiseExp(blockNumber, pars, mess)

    if ismember(blockNumber, [0, 1])
        session.runs       = generateRuns(pars);
        session.currentN   = 1;
        session.subject    = input(sprintf('[%s: ', mess.subj), 's'); 
        fprintf([char(8), ']\n']);
        session.sessN      = str2num(input(sprintf('[%s: ', mess.sessn), 's'));
        fprintf([char(8), ']\n']);
        save(pars.cache, '-struct', 'session');
    elseif blockNumber < 0
        fprintf('\nResuming from cache! ')
        session = load(pars.cache);
        fprintf('Subject "%s", ', session.subject);
        fprintf('session #%d. ', session.sessN)
        fprintf('Last block cached: block #%d\n\n', session.currentN - 1);
    elseif mod(blockNumber, 1) == 0 && blockNumber <= length(pars.scheme)
        fprintf('Attemting to read block designs from cache...\n')
        try 
            session = load(pars.cache);
            if session.currentN ~= blockNumber
                warning('Cached block number differs from specified number!')
                fprintf('Last block cached: block #%d\n\n', session.currentN-1);
            end
        catch
            warning('No cache found! Generating new block desing!\n')
            session.subject = input(sprintf('[%s: ', mess.subj), 's');
            fprintf([char(8), ']\n']);
            session.sessN      = str2num(input(sprintf('[%s: ', mess.sessn), 's'));
            fprintf([char(8), ']\n']);
            session.runs = generateRuns(pars);
        end
        session.currentN   = blockNumber;
        save(pars.cache, '-struct', 'session');
    else
        error('Resume option not valid! (use <0 or integers <= %d)!', ...
        												length(pars.schema));
    end

    estimateTime(pars);

end



function estimateTime(pars)

    nOfOddBallRuns = 0; nOfTonotopRuns = 0; nOfFuncLocRuns = 0;
    for i = pars.scheme
        nOfOddBallRuns = nOfOddBallRuns + strcmp(pars.task{i}, 'omissions');
        nOfFuncLocRuns = nOfFuncLocRuns + strcmp(pars.task{i}, 'localiser');
    end

    % 1. OMISSION RUNS

    % 1.1 Time per trial
    tTone   = sum(pars.ob.dur) + pars.ob.ISI;
    tTrain  = mean(pars.ob.len) * tTone - pars.ob.ISI;  

    % 1.3 Expected ITI (from a truncated distribution)
    alph2  = (pars.ob.feedbackT + tTrain - pars.ob.avgEOA) / sqrt(2);
    beth2  = (pars.ob.maxEOA) / sqrt(2);
    muDelt = (exp(-alph2 ^ 2) - exp(-beth2 ^ 2)) / (erf(beth2) - erf(alph2));
    avgEOA = pars.ob.avgEOA + sqrt(2/pi) * muDelt;
    avgITI = avgEOA - tTrain;

    % 1.4 Time on oddball trials
    trTrials = length(pars.ob.devLocs) * length(pars.ob.tones);
    trTrials = trTrials * (1 + pars.ob.probNoOm); % + non-omission trials
    tmTrials = trTrials * (tTrain + avgITI + 0.5 * pars.ob.feedbackT) - avgITI;

    % 1.5 Time on null events
    nOfNullEvs = length(pars.ob.tones)  * (1 + pars.ob.probNoOm) * pars.ob.null;
    avgITI   = pars.ob.avgEOA - sqrt(2/pi) * exp(- alph2 ^ 2) / (1 - erf(alph2));
    nullTime = (tTrain + avgITI) * nOfNullEvs;

    % 1.6 Time per run 
    tOddBallRun = tmTrials + nullTime;


    % 2. FUNCLOC RUNS

    % 2.1 Expected ITI (from a truncated distribution)
    alph = (pars.loc.minIBI - pars.loc.avgIBI);
    beth = (pars.loc.maxIBI - pars.loc.avgIBI);
    A = (1 / sqrt(2*pi)) * (exp( -0.5 * alph^2) - exp( -0.5 * beth^2));
    Z = (1 / 2)          * (erf(beth / sqrt(2)) - erf(alph / sqrt(2)));
    expectedIBI = pars.loc.avgIBI + A / Z;
    if isnan(expectedIBI) % This happens when max = avg = min (dist not defined)
        expectedIBI = pars.loc.avgIBI
    end

    % 2.2 Expected time per funcloc block
    timePerBlock = (pars.loc.dur + pars.loc.ISI) * pars.loc.blockLen + expectedIBI;
    tFuncLocRun = (timePerBlock * 2 * pars.loc.Nblocks - expectedIBI);

    % 3. TOTAL TIME
    totalTime = nOfOddBallRuns * tOddBallRun + ...
                nOfFuncLocRuns * tFuncLocRun; 
    totalTime = totalTime + (pars.IRI + 20) * (length(pars.scheme) - 1);

    pr = [sprintf('%s -> ', pars.task{pars.scheme}), repmat(char(8), [1, 4])];
    fprintf('[Programme: %s]\n', pr)
    fprintf('[Estimated times: oddBallRun = %.1f;', tOddBallRun / 60); 
    fprintf(' funcLocRun = %.1fm;', tFuncLocRun / 60);
    fprintf(' experiment = %.1fm]\n', totalTime / 60);

end



function handles = initialiseHandles(pars, subject)

    [machine, fMRI] = checkMachine();

    if not(fMRI), fprintf('[WARNING: fmri not detected]\n'); end
    if fMRI, fprintf('[fmri detected]\n'); end

    addPsychToolBoxPath();

    % Checking everything is loaded and ready
    KbCheck;
    WaitSecs(0.1);
    GetSecs;
    InitializePsychSound(1);

    Screen('CloseAll');      % Just in case
    PsychPortAudio('Close'); %

    audioPort = openAudioPort(pars);

    verboseFlagHandle = Screen('Preference', 'VisualDebuglevel', 3);

    % Skips the sync testing, which can take up to 3s (good for testing)
    Screen('Preference','SkipSyncTests', 1);

    screenID     = 0; % Presentation screen usually last 
    bckgCol      = [0 0 0];                % Black RGB
    if fMRI
        screenLoc    = [480 -600 1280 0] ;  % Pixel coordinates for the window
    else
        screenLoc    = [0 0 800 600];
    end
    windowH = Screen('OpenWindow', screenID, bckgCol, screenLoc);
    
    %if fMRI
    %     windowH = Screen('OpenWindow', screenID, bckgCol);
    %else
         %windowH = Screen('OpenWindow', screenID, bckgCol, screenLoc);
    %end
    
    Screen('TextFont', windowH, 'Arial');   
    Screen('TextColor', windowH, [0 0 0]); % Sets the text to RGB [,,]
    Screen('TextSize', windowH, 40);
    % Screen(windowH,'Flip');
    % HideCursor;

    handles.window = windowH;
    handles.verb   = verboseFlagHandle;
    handles.audio  = audioPort;

    % Objects keys/keysTrig mappings key labels to actual keys / triggers
    handles.keys.trigger   = KbName('t'); % Trio scanning trigger 
    handles.keys.space     = KbName('space');
    handles.keys.button(1) = KbName('z');
    handles.keys.button(2) = KbName('g');
    handles.keys.button(3) = KbName('r');
    handles.keys.button(4) = KbName('m');
    if isunix
        handles.keys.escape = KbName('escape');
    else
        handles.keys.escape = KbName('esc');
    end
    
end



function [logFile, messFile] = openLogfiles(session, n, pars)

    if exist('./logFiles') == 0
        mkdir('logFiles')
    end

    prefix = sprintf('logFiles/sub%s-sess%02d',session.subject,session.sessN);
    task   = pars.task{pars.scheme(n)};

    % Opening logFile
    logFileName = sprintf('%s-run%02d-%s-events.tsv', prefix, n, task);
    if exist(logFileName, 'file') == 2
        warning('Log file %s already existed and was renamed!', logFileName)
        movefile(logFileName, sprintf('%s_old%d', logFileName, randi(99)));
    end
    logFile = fopen(logFileName, 'a');

    messFileName = sprintf('%s-run%02d-%s-mess.tsv', prefix, n, task);
    if exist(messFileName, 'file') == 2
        warning('Mess file %s already existed and was renamed!', messFileName)
        movefile(messFileName, sprintf('%s_old%d', messFileName, randi(99)));
    end
    messFile = fopen(messFileName, 'a');

    stamp = string(datetime('now','TimeZone','local', ...
                                  'Format','d-MMM-y HH:mm:ss'));
    sessInf = sprintf('Sub: %s, Sess: %d, Run %d', ...
                                  session.subject, session.sessN, n);

    meta = sprintf('%s logfile. %s. Creation date: %s\n', task,sessInf,stamp);
    fprintf(logFile, meta);

    meta = sprintf('%s messfile. %s. Creation date: %s\n', task,sessInf,stamp);
    fprintf(messFile, meta);

    if strcmp(task, 'omissions')
        logHead = sprintf('onset\tduration\ttype\tanswer\tRT\n');
    elseif strcmp(task, 'localiser')
        logHead = sprintf('onset\tduration\ttype\tcorrect\n');
    end
    fprintf(logFile, logHead);
    
    messHead = sprintf('onset\tmessage\n');
    fprintf(messFile, messHead);

end



function triggerTime = waitForScanner(handles, mess)

    displayText(handles, mess.waiting, false);

    keyPressTime = [];

    [~, triggerTime, keyCode] = KbCheck();
    
    while (keyCode(handles.keys.trigger) == 0)
        [~, triggerTime, keyCode] = KbCheck(); 
        if (keyCode(handles.keys.escape))
            terminateExecution();
        end
    end

    consciousWait(0.5, handles);

    displayText(handles, '+', false);
    
    for i = 1:4, % Wait for 4 scans before starting the experiment
        [~, ~, keyCode] = KbCheck();
        while (keyCode(handles.keys.trigger) == 0)
            [~, ~, keyCode] = KbCheck(); 
            if (keyCode(handles.keys.escape))
                terminateExecution();
            end
        end
        consciousWait(0.5, handles); 
    end

end



function waitForSpace(handles)

    [~, ~, keyCode] = KbCheck();
    
    while (keyCode(handles.keys.space) == 0)
        [~, ~, keyCode] = KbCheck(); 
        if (keyCode(handles.keys.escape))
            terminateExecution(handles);
        end
    end

end



function [keyPressed, keyPressTime] = waitForKey(handles, timeOut)

    keyPressTime = [];
    t0 = GetSecs();
    
    [~, keyPressTime, keyCode] = KbCheck();
    
    while (sum(keyCode(handles.keys.button))==0) && (GetSecs()-t0 < timeOut)
        [~, keyPressTime, keyCode] = KbCheck(); 
        if (keyCode(handles.keys.escape))
            terminateExecution(handles);
        end
    end

    if isempty(keyPressTime)
        keyPressed = [];
    else   
        KbReleaseWait();
        keyPressed = KbName(keyCode);
    end

end



function button = keyToButton(keyPressed, handles)

    button = -1;
    keys = ['z' 'g' 'r' 'm'];
    for i = 1:length(handles.keys.button)
        if keyPressed == keys(i)
            button = i;
        end
    end

end



function consciousWait(waitDuration, handles)

    if nargin == 1
        if isunix
            escape = KbName('escape');
        else
            escape = KbName('esc');
        end
    else
        escape = handles.keys.escape;
    end

    t0 = GetSecs();

    while (GetSecs() - t0 < waitDuration)
        [~, ~, keyCode] = KbCheck(); 
        if keyCode(escape)
            terminateExecution();     
        end
    end

end



function writeTrainToLog(handles, onsets, blck, comb, answer, rt, pars)

    if blck(2) > -1
        dur    = sum(pars.ob.dur);
        trID   = sprintf('loc%d-%s', blck(2), comb);
 
        logInput = [];
        for i = 1:length(onsets)
            if i == 1
                trialType = sprintf('%s-st0', trID);
            elseif i == blck(2)
                trialType = sprintf('%s-om', trID);
            else
                trialType = sprintf('%s-std', trID);
            end
            logLine = sprintf('%.4f\t%.4f\t%s', onsets(i), dur, trialType);
            if i == blck(2)
                if rt == []
                    logLine = sprintf('%s\t%d\t NaN', logLine, answer);
                else
                    logLine = sprintf('%s\t%d\t%.4f', logLine, answer, rt);
                end
            end
            logInput = [logInput, logLine, '\n'];
        end
    else
        dur = (sum(pars.ob.dur) + pars.ob.ISI) * pars.ob.len - pars.ob.ISI;
        logInput = sprintf('%.4f\t%.4f\tnull\n', onsets(1), dur);
    end
    
    fprintf(handles.logF, logInput);

end 



function writeSoundToLog(handles, onset, dur, tKey, result, trialID)

    if strcmp(result, 'falseAlarm') || strcmp(result, 'miss')
        correct = -1; % error
    elseif strcmp(result, 'detection')
        correct = 1;  % detection
    else
        correct = 0;  % true positive
    end

    L = sprintf('%.4f\t%.4f\t%s\t%d\n', onset, dur, trialID, correct);
    fprintf(handles.logF, L);

    if (tKey > 0)
        L = sprintf('%.4f\t%.4f\t%s\t%s\n', tKey, 0.1, 'key', result);
        fprintf(handles.logF, L);
    end

end



function extraWait = computeITI(blck, i, pars, rt)

    if or(rt ==-1, isempty(rt)), rt = pars.ob.maxRT; end
    if rt < 0, rt = 0; end
    if blck(i, 2) == 0,
        alreadyWaited = (pars.ob.len-6)*(sum(pars.ob.dur)+pars.ob.ISI);
        ITI = pars.ob.minITI + abs(0.1*randn());
    else 
        timeSincePrev = (pars.ob.len-blck(i,2))*(sum(pars.ob.dur)+pars.ob.ISI);
        alreadyWaited = max(0, rt + pars.ob.feedbackT - timeSincePrev);
        nextLoc = blck(min(i + 1, size(blck, 1)), 2);
        timeToNext = (nextLoc - 1) * (pars.ob.ISI + sum(pars.ob.dur));
        EOA = min(pars.ob.maxEOA, max(pars.ob.minEOA, pars.ob.avgEOA + randn()));
        ITI = max(pars.ob.minITI, EOA - timeSincePrev - timeToNext);
    end

    extraWait = max(0, ITI - alreadyWaited);

end



function displayText(handles, message, saveToLog)

    if nargin < 3, saveToLog = true; end;

    if ~handles.window == 0
        if strcmp(message, '+red')
            textColor = [200, 65, 0];
            mess = '+';
        elseif strcmp(message, '+green') 
            textColor = [65, 200, 0];
            mess = '+';
        else
            textColor = [255, 255, 255];
            mess = message;  
        end
        DrawFormattedText(handles.window, mess, 'center', 'center', textColor);
        Screen(handles.window, 'Flip');
    end

    if saveToLog
        mess = sprintf('%.4f\t%s\n', GetSecs() - handles.absOnset, message);
        fprintf(handles.messF, mess);
    end

    if not(strcmp(message, '+'))
        fprintf(['--Screen: ', message, '\n']);
    end

end



function soundCheck(pars)

    addPsychToolBoxPath();

    % Checking everything is loaded and ready
    InitializePsychSound(1);
    audioPort = openAudioPort(pars);
    s = zeros([1, 0.05 * pars.fs]);

    if ismember('omissions', pars.task(pars.scheme))
        for j = 1:5:length(pars.ob.tones)
            stim = createSinusoid(pars.ob.tones(j), pars);
            s = [s, stim.s, zeros([1, 0.2 * pars.fs])];
        end
    end

    s = [s, zeros([1, 0.5*pars.fs])];

    if ismember('localiser', pars.task(pars.scheme))
        stim = [];
        sounds = readSounds(pars);
        for ix = [1 3 14]
            stim  = [stim, pars.natsoundGain * sounds{ix}.waveform];
        end
        s = [s, stim, zeros([1, 0.2*pars.fs])];
    end

    PsychPortAudio('FillBuffer', audioPort, s(1:(length(s) - pars.fs)));
    PsychPortAudio('Start', audioPort); 
    PsychPortAudio('Stop', audioPort, 1, 1);

end



function findPort(pars)

    fprintf('\nTesting all available soundports...\n\n');
    addPsychToolBoxPath();

    % Checking everything is loaded and ready
    InitializePsychSound(1);
    PsychPortAudio('Close'); %

    stim = [];
    sounds = readSounds(pars);
    s = [pars.natsoundGain * sounds{1}.waveform,  zeros([1, pars.fs])];

    devs = PsychPortAudio('GetDevices');

    fprintf('\n')
    for i = 1:length(devs)
        fprintf('devIx %d -> devName "%s" ...', devs(i).DeviceIndex, devs(i).DeviceName);
        try
            audioPort = openAudioPort(pars, devs(i).DeviceIndex);
            PsychPortAudio('FillBuffer', audioPort, s(1:(length(s) - pars.fs)));
            PsychPortAudio('Start', audioPort); 
            PsychPortAudio('Stop', audioPort, 1, 1);
            consciousWait(1);
        catch exception
            warning('Error was thrown while trying to use this port...');
            fprintf(getReport(exception));
        end
        fprintf('\n');
    end

    closeInterfaces();

end



function audioPort = openAudioPort(pars, devIx)

    devs = PsychPortAudio('GetDevices');
    if nargin == 1
        
        [machine, fMRI] = checkMachine();
    
        if fMRI
            devs = PsychPortAudio('GetDevices');
            devName = 'Microsoft Sound Mapper - Output'; 
            for i = 1:length(devs)
                if strcmp(devs(i).DeviceName, devName)
                    devIx = devs(i).DeviceIndex;
                end
            end
        elseif strcmp('ornette', machine)
            devName = 'HDA Intel PCH: ALC887-VD Analog (hw:1,0)'; 
            for i = 1:length(devs)
                if strcmp(devs(i).DeviceName, devName)
                    devIx = devs(i).DeviceIndex;
                end
            end
        else
            devIx = [];
        end
    end

    PsychPortAudio('Close');
    ctrlc = 'com.mathworks.mde.cmdwin.CmdWinMLIF.getInstance().processKeyFromC(2,67,''C'')';
    killTimer = timer('TimerFcn', ctrlc, 'StartDelay', 3);
    start(killTimer);
    audioPort = PsychPortAudio('Open', devIx, 1, 2, pars.fs, 1, 0);
    stop(killTimer);

end



function addPsychToolBoxPath()

    [machine, fMRI] = checkMachine();

    if exist('PsychPortAudio') == 0
        if strcmp('ornette', machine) || strcmp('intrepid', machine)
            addpath(genpath('/usr/share/psychtoolbox-3'));
        elseif strcmp('coldplay', machine)
            softPath = '/afs/cbs.mpg.de/software/psychtoolbox/';
            psyVers  = '0.20170103/ubuntu-xenial-amd64/9.1/Psychtoolbox/';
            addpath(genpath([softPath, psyVers])); 
        elseif fMRI
            softPath = 'C:\Program Files\toolbox\Psychtoolbox';
            psyVers  = '0.20170103/ubuntu-xenial-amd64/9.1/Psychtoolbox/';
        else
            error('Psychtoolbox path not specified for %s!', machine)
        end
    end

end



function [machine, fMRI] = checkMachine()

    [~, machine] = system('hostname');
    machine = string(machine(1:(length(machine)-1)));
    fMRI = strcmp('pc10present3', machine);   
    
end
    


function closeInterfaces()

    PsychPortAudio('Close');
    Screen('CloseAll');
    ShowCursor;
    fclose('all');

end



function terminateExecution()

    closeInterfaces();
    error('Execution terminated by user (not a real error!)');

end



function sounds = readSounds(pars)

    % 1. Gather all the sounds in the stimuli path
    soundDir  = dir([pars.loc.stimPath '/*s3_*.wav']);
    
    % Now we process each of the sounds independently
    for i = 1:length(soundDir)
        % 2. We write the wavefile to x and the sampling rate to fs
        [x, fs] = audioread([pars.loc.stimPath '/' soundDir(i).name]);

        % 3. If the sounds' (fs) and our desired sampling rate (pars.fs) 
        %    differ, we resample the sound to pars.fs (don't worry about this)
        if (fs ~= pars.fs)
            if (fs == 16000) && (pars.fs == 48000)
                q = 1; p = 3;
            elseif (fs == 44100) && (pars.fs == 48000)
                q = 147; p = 160;
            else
                q = fs / 100; p = pars.fs / 100;
            end
            x = resample(x, p, q);
        end

        % 4. Prepare the label for the sound based ont he filename
        name = strsplit(soundDir(i).name, '_');

        % 5. We store the (resampled) sound's waveform in the "sounds" array
        sounds{i}.waveform = x';
        sounds{i}.label    = [name{2:3}];
    end

end



function mess = getMessages(pars)

    mess.pressEnter = 'Press ENTER to continue';
    mess.trial      = 'Seq';
    mess.waiting    = 'Warten auf den Start...';
    mess.waitRun    = 'PAUSE';
    mess.endExp     = 'Fertig! Vielen dank! :)';
    mess.subj       = 'Subject reference';
    mess.sessn      = 'Session number';

end



function pars = getParameters()
    
    % Global parameters
    pars.cache   = 'sessionCache.mat';
    pars.fs      = 48000; % Herz
    pars.IRI     = 60;  % Inter-run-interval, seconds
    pars.scheme  = [1 1 1 2 1];
    pars.task    = {'omissions', 'localiser'};

    % Oddball parameters
    ob.dur       = [0.005, 0.04, 0.005];
    ob.devLocs   = [4:6]; % possible locations of the deviants
    ob.tones     = logspace(log10(500), log10(1500), 21);
    ob.len       = 8;   % train length
    ob.ISI       = 0.5; % seconds
    ob.maxEOA    = 11;  % seconds
    ob.minEOA    = 3;   % seconds
    ob.minITI    = 2;   % seconds
    ob.avgEOA    = 5;   % seconds 
    ob.maxRT     = 2;   % seconds
    ob.feedbackT = 1.5; % seconds
    ob.null      = 0.3; % additive proportion of null events (not total!)
    ob.probNoOm  = 2/(4*3); % ~0.167; Must be = integer / (nOfRuns * devlocs)
    ob.standards = 1:length(ob.tones); 

    % Localiser parameters
    loc.dur      = 1;    % second
    loc.ISI      = 0.1;  % seconds
    loc.minIBI   = 0;    % seconds
    loc.avgIBI   = 0;    % seconds
    loc.maxIBI   = 0;    % seconds
    loc.maxRT    = loc.dur + loc.ISI; % seconds
    loc.repProb  = 0.05; % The probability of a sound being repeated
    loc.blockLen = 16;   % Number of sounds in each block
    loc.Nblocks  = 10;   % Number of silence/sound blocks
    loc.stimPath = './stimuli/';

    pars.phons        = 80;
    globalGain        = 0.8 * 1/loudWeight(min(ob.tones), pars.phons);
    pars.beepGain     = 0.4 * globalGain;
    pars.natsoundGain = globalGain;

    pars.ob  = ob;
    pars.loc = loc;

end
