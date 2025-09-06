%% JANUARY 2025

clear all
close all
clc

%% User inputs 
player = 0;                             % You! 0: white; 1: black
depth = 2;                              % Analysis depth (moves per side)
max_Time = 1800;                        % Max thinking seconds per move
TwoEngines = 0;                         % 0: you play; 1: engine vs engine
sideToShowWhenTwoEngines = 0;           % Show white/black side (0/1)
A0 = [14, 13, 12, 15, 16, 12, 13, 14;   % Initial board configuration
      11, 11, 11, 11, 11, 11, 11, 11;
      00, 00, 00, 00, 00, 00, 00, 00;
      00, 00, 00, 00, 00, 00, 00, 00;
      00, 00, 00, 00, 00, 00, 00, 00;
      00, 00, 00, 00, 00, 00, 00, 00;
      01, 01, 01, 01, 01, 01, 01, 01;
      04, 03, 02, 05, 06, 02, 03, 04];

% 00 = empty
% 01 = White Pawn
% 02 = White Bishop
% 03 = White kNight
% 04 = White Rook
% 05 = White Queen
% 06 = White King
% 11 = Black Pawn
% 12 = Black Bishop
% 13 = Black kNight
% 14 = Black Rook
% 15 = Black Queen
% 16 = Black King

MfileName = "ex29.mat";                  % Output file with the game
startfromthisconfiguration = 0;          % 1 to continue an old game
ignoreLastNMoves = 0;                    % N of last moves to ignore
thisconfiguration = "ex29.mat";          % Old game to load (here as Black)

%% Utility and internal variables
if(TwoEngines==1)
    player = 0;
end
global max_iter
global max_time
global t1
global t2
global mvv
global prev
global sside
sside = sideToShowWhenTwoEngines;
mvv = 0;
prev = 0;
max_iter = depth;
max_time = max_Time;
max_time = max(5.0, (max_time/2.0-5.0));
t1 = 0;
t2 = 0;
global whitePiecesColor
global blackPiecesColor
global whiteBackgroundColor
global blackBackgroundColor
whitePiecesColor = [1,1,1];
blackPiecesColor = [0,0,0];
whiteBackgroundColor = [0.8,0.8,0.5];
blackBackgroundColor = [0.7,0.4,0.1];
global playAlone;
playAlone = TwoEngines;

% Struct containing all the information about a certain move
M(1).A = A0;
M(1).turn = 0;                     % 0 = white; 1 = black.
M(1).whiteAlreadyCastled = 0;
M(1).whiteLeftRookMoved = 0;
M(1).whiteRightRookMoved = 0;
M(1).whiteKingMoved = 0;
M(1).blackAlreadyCastled = 0;
M(1).blackLeftRookMoved = 0;
M(1).blackRightRookMoved = 0;
M(1).blackKingMoved = 0;
M(1).lastMove = [0,0,0,0];         % x1, y1, x2, y2
checkStatus0.out_check = 0;
checkStatus0.out_zi = 0;
checkStatus0.out_zj = 0;
M(1).checkStatus = checkStatus0;
M(1).checkmateStatus = 0;
M(1).staleStatus = 0;
M(1).captureOrPawnMove = 0;
M(1).staleRepetitions = 0;
M(1).staleMoves = 0;
M(1).player = player;
M(1).Depthch = 0;
M(1).thinkingTime = 0;

if(startfromthisconfiguration==1)
    load(thisconfiguration);
    M = M(1:(end-ignoreLastNMoves));
end

%% Play the game
plotConfig(M(end));
if(player==0)
    ply="WHITE";
else
    ply="BLACK";
end
text(520,-250, strcat("Welcome! Click anywhere to start playing as ",...
    ply, "."), 'FontSize',12)
waitforbuttonpress;
prev = 1;

while(1)
    % White turn
    A = M(end).A;
    M(end).checkStatus = isCheck(M(end));
    M(end).checkmateStatus = isCheckmate(M(end));
    M(end).staleStatus = isStale(M(end));
    M(end).staleRepetitions = isStaleRepetitions(M);
    M(end).staleMoves = isStaleMoves(M);
    if(M(end).checkmateStatus == 1 || M(end).staleStatus==1 ||...
            M(end).staleRepetitions==1 || M(end).staleMoves==1)
        break;
    end
    if(playAlone == 1)
        player=1-player;
        M(end).player = player;
    end
    if(player==0)
        S_out = listNewPossibleConfigurations(M(end));
        isValidMove = 0;
        t9 = tic();
        while(isValidMove==0)
            try
                plotConfig(M(end));
                pause(1.0);
                [i1, j1] = collectUserMove(M(end));
                [i2, j2] = collectUserMove(M(end));
                for ii=1:length(S_out)
                    if(sum(S_out(ii).M.lastMove == [i1,j1,i2,j2])==4)
                        isValidMove = 1;
                        haveToChoice = 0;
                        if(ii<length(S_out))
                            for jj=(ii+1):length(S_out)
                                if(sum(S_out(jj).M.lastMove == ...
                                        [i1,j1,i2,j2])==4)
                                    haveToChoice = 1;
                                    break;
                                end
                            end
                        end
                        if(haveToChoice == 0)
                            M(end+1) = S_out(ii).M;
                            M(end).Depthch = 0;
                            M(end).thinkingTime = toc(t9);
                        else
                            [uy, validChoice] = askChoice(M(end).turn);
                            isValidMove = validChoice;
                            if(isValidMove==1)
                                M(end+1) = S_out(ii+uy-2).M;
                                M(end).Depthch = 0;
                                M(end).thinkingTime = toc(t9);
                            end
                        end
                        break;
                    end
                end
            catch
                isValidMove=0;
            end
            if(isValidMove==0)
                title("NOT A VALID MOVE: TRY AGAIN!")
                pause(1);
                title("");
            end
        end
    else
        plotConfig(M(end));
        pause(1.0);
        t5=tic();
        [Mch, Depthch] = chooseBestMove(M, depth);
        M(end+1) = Mch;
        t6=toc(t5);
        M(end).Depthch = Depthch;
        M(end).thinkingTime = t6;
        int32([Depthch, t6])
        pause(1);
    end
    mvv = mvv + 1;

    % Black turn
    A = M(end).A;
    M(end).checkStatus = isCheck(M(end));
    M(end).checkmateStatus = isCheckmate(M(end));
    M(end).staleStatus = isStale(M(end));
    M(end).staleRepetitions = isStaleRepetitions(M);
    M(end).staleMoves = isStaleMoves(M);
    if(M(end).checkmateStatus == 1 || M(end).staleStatus==1 || ...
            M(end).staleRepetitions==1 || M(end).staleMoves==1)
        break;
    end
    if(playAlone == 1)
        player=1-player;
        M(end).player = player;
    end
    if(player==0)
        plotConfig(M(end));
        pause(1.0);
        t5=tic();
        [Mch, Depthch] = chooseBestMove(M, depth);
        M(end+1) = Mch;
        t6=toc(t5);
        M(end).Depthch = Depthch;
        M(end).thinkingTime = t6;
        int32([Depthch, t6])
        pause(1);
    else
        S_out = listNewPossibleConfigurations(M(end));
        isValidMove = 0;
        t9 = tic();
        while(isValidMove==0)
            try
                plotConfig(M(end));
                pause(1.0);
                [i1, j1] = collectUserMove(M(end));
                [i2, j2] = collectUserMove(M(end));
                for ii=1:length(S_out)
                    if(sum(S_out(ii).M.lastMove == [i1,j1,i2,j2])==4)
                        isValidMove = 1;
                        haveToChoice = 0;
                        if(ii<length(S_out))
                            for jj=(ii+1):length(S_out)
                                if(sum(S_out(jj).M.lastMove == ...
                                        [i1,j1,i2,j2])==4)
                                    haveToChoice = 1;
                                    break;
                                end
                            end
                        end
                        if(haveToChoice == 0)
                            M(end+1) = S_out(ii).M;
                            M(end).Depthch = 0;
                            M(end).thinkingTime = toc(t9);
                        else
                            [uy, validChoice] = askChoice(M(end).turn);
                            isValidMove = validChoice;
                            if(isValidMove==1)
                                M(end+1) = S_out(ii+uy-2).M;
                                M(end).Depthch = 0;
                                M(end).thinkingTime = toc(t9);
                            end
                        end
                        break;
                    end
                end
            catch
                isValidMove = 0;
            end
            if(isValidMove==0)
                title("NOT A VALID MOVE: TRY AGAIN!")
                pause(1);
                title("");
            end
        end
    end
end

A = M(end).A;
plotConfig(M(end));
save(MfileName, "M");

%% Functions

% Main recursive algorithm
function [bestMove, pts, why] = findBestMove(M_vec, it)
    global max_iter
    global max_time
    global t1
    global t2
    if(it<0.5)
        t1 = tic;
    end
    pts = 0;
    M = M_vec(end);
    why = 0;

    % Exit case: last iteration
    if(it>0.5)
        if(isStaleRepetitions(M_vec)==1)
            bestMove = M;
            pts = 0;
            return;
        elseif(isStaleMoves(M_vec)==1)
            bestMove = M;
            pts = 0;
            return;
        elseif(isCheckmate(M)==1)
            bestMove = M;
            if(M.player==0)
                pts = +1000000;
            else
                pts = -1000000;
            end
            return;
        elseif(isStale(M)==1)
            bestMove = M;
            pts = 0;
            return;
        elseif(it>=max_iter)
            bestMove = M;
            pts = pointsOfM(M);
            return;
        end
    end
    
    % Generic iteration
    S_outA = listNewPossibleConfigurations(M);
    L_S_outA = length(S_outA);
    %S_outA = S_outA(randperm(L_S_outA));
    bestMove = S_outA(1).M;
    pts = pointsOfM(S_outA(1).M);
    bestPossible=0;
    
    for ii=1:L_S_outA
        t2 = toc(t1);
        if(t2>max_time)
            why=1;
            break;
        end
        Ma = S_outA(ii).M;
        S_outB = listNewPossibleConfigurations(Ma);
        L_S_outB = length(S_outB);
        if(L_S_outB==0)
            M_vec_provv = [M_vec,Ma];
            if(isStaleRepetitions(M_vec_provv)==1)
                Ypts = 0;
            elseif(isStaleMoves(M_vec_provv)==1)
                Ypts = 0;
            elseif(isCheckmate(Ma)==1)
                if(M.player==0)
                    Ypts = -500000+it/50;
                else
                    Ypts = +500000-it/50;
                end
            elseif(isStale(Ma)==1)
                Ypts = 0;
            end
            Sbests(ii).best = Ma;
            Sbests(ii).pts = Ypts;
            vecpts(ii) = Ypts;
            continue;
        end
        
        %S_outB = S_outB(randperm(L_S_outB));
        if(M.player==0)
            worstpts = -1000000;
        else
            worstpts = +1000000;
        end
        for jj=1:L_S_outB
            Mb = S_outB(jj).M;
            M_vec_copy = [M_vec, Ma, Mb];
            [~, pts_jj, ~] = findBestMove(M_vec_copy, it+1);
            if(M.player==0)
                if(pts_jj>worstpts)
                    worstpts = pts_jj;
                end
            else
                if(pts_jj<worstpts)
                    worstpts = pts_jj;
                end
            end
        end

        Sbests(ii).best = Ma;
        if(M.player==0)
            Sbests(ii).pts = worstpts+it/50;
            vecpts(ii) = worstpts+it/50;
        else
            Sbests(ii).pts = worstpts-it/50;
            vecpts(ii) = worstpts-it/50;
        end
        
    end

    gg = 0;
    if(ii>1.5 || bestPossible==1)
        dmm = 0.01;
        if(M.player==0)
            minspts = min(vecpts);
            for uu=1:length(Sbests)
                if(Sbests(uu).pts<(minspts+dmm))
                    gg = gg+1;
                    possible_pts(gg) = Sbests(uu).pts;
                    possible_bestMove(gg) = Sbests(uu).best;
                end
            end
        else
            maxspts = max(vecpts);
            for uu=1:length(Sbests)
                if(Sbests(uu).pts>(maxspts-dmm))
                    gg = gg+1;
                    possible_pts(gg) = Sbests(uu).pts;
                    possible_bestMove(gg) = Sbests(uu).best;
                end
            end
        end
    end
    if(gg>0.5)
        choicee = randi(gg);
        pts = possible_pts(choicee);
        bestMove = possible_bestMove(choicee);
    end
end

% Evaluate the score of a given board configuration
function [value, pointsWhite, pointsBlack] = pointsOfM(M)
    pointsWhite = 0;
    pointsBlack = 0;
    A = M.A;
    bgBlack = [0,0];
    bgWhite = [0,0];
    iWk = 0;
    jWk = 0;
    iBk = 0;
    jBk = 0;
    dWhiteToBlackKing = 0;
    dBlackToWhiteKing = 0;
    for i=1:8
        for j=1:8
            if(A(i,j)==1)
                pointsWhite = pointsWhite + 1;
            elseif(A(i,j)==2)
                pointsWhite = pointsWhite + 3;
            elseif(A(i,j)==3)
                pointsWhite = pointsWhite + 3;
            elseif(A(i,j)==4)
                pointsWhite = pointsWhite + 5;
            elseif(A(i,j)==5)
                pointsWhite = pointsWhite + 9;
            elseif(A(i,j)==6)
                pointsWhite = pointsWhite + 1000;
            elseif(A(i,j)==11)
                pointsBlack = pointsBlack + 1;
            elseif(A(i,j)==12)
                pointsBlack = pointsBlack + 3;
            elseif(A(i,j)==13)
                pointsBlack = pointsBlack + 3;
            elseif(A(i,j)==14)
                pointsBlack = pointsBlack + 5;
            elseif(A(i,j)==15)
                pointsBlack = pointsBlack + 9;
            elseif(A(i,j)==16)
                pointsBlack = pointsBlack + 1000;
            end
            bgBlack = bgBlack + [i,j]*(A(i,j)>10);
            bgWhite = bgWhite + [i,j]*(A(i,j)<10 && A(i,j)>0.5);
            if(A(i,j)==6)
                iWk = i;
                jWk = j;
            elseif(A(i,j)==16)
                iBk = i;
                jBk = j;
            end
        end
    end
    bgBlack = bgBlack/sum(A>10,"all");
    bgWhite = bgWhite/sum((A<10&A>0.5),"all");
    if(iWk>0.5&&jWk>0.5)
        dBlackToWhiteKing = sqrt((bgBlack(1)-iWk)^2 + (bgBlack(2)-jWk)^2);
    else
        dBlackToWhiteKing = 0;
    end
    if(iBk>0.5&&jBk>0.5)
        dWhiteToBlackKing = sqrt((bgWhite(1)-iBk)^2 + (bgWhite(2)-jBk)^2);
    else
        dWhiteToBlackKing = 0;
    end
    dBlackToCenter = sqrt((bgBlack(1)-4.5)^2 + (bgBlack(2)-4.5)^2);
    dWhiteToCenter = sqrt((bgWhite(1)-4.5)^2 + (bgWhite(2)-4.5)^2);
    dWhiteToPromotion = abs(bgWhite(1)-1);
    dBlackToPromotion = abs(bgBlack(1)-8);
    if(pointsWhite<=1002)
        pointsWhite = pointsWhite - dWhiteToPromotion;
    end
    if(pointsBlack<=1002)
        pointsBlack = pointsBlack - dBlackToPromotion;
    end
    if(M.whiteAlreadyCastled==1)
        pointsWhite = pointsWhite + 1;
    end
    if(M.blackAlreadyCastled==1)
        pointsBlack = pointsBlack + 1;
    end
    if(M.whiteKingMoved==1)
        pointsWhite = pointsWhite - 0.3;
    end
    if(M.blackKingMoved==1)
        pointsBlack = pointsBlack - 0.3;
    end
    if(M.whiteLeftRookMoved==1)
        pointsWhite = pointsWhite - 0.02;
    end
    if(M.whiteRightRookMoved==1)
        pointsWhite = pointsWhite - 0.02;
    end
    if(M.blackLeftRookMoved==1)
        pointsBlack = pointsBlack - 0.02;
    end
    if(M.blackRightRookMoved==1)
        pointsBlack = pointsBlack - 0.02;
    end

    pointsWhite = pointsWhite - dWhiteToBlackKing/10.0 - ...
        dWhiteToCenter/100.0;
    pointsBlack = pointsBlack - dBlackToWhiteKing/10.0 - ...
        dBlackToCenter/100.0;

    value = pointsWhite - pointsBlack;

end

% Evaluate whether it's stalemate for no available moves
function [staleStatus] = isStale(M)
    staleStatus = 0;
    checkStatus = isCheck(M);
    out_check = checkStatus.out_check;
    if(out_check==1)
        return;
    end
    S_out = listNewPossibleConfigurations(M);
    if(isempty(S_out))
        staleStatus = 1;
    end

end

% Evaluate whether it's stalemate for 3 repetitions of the same board
function [outRepet] = isStaleRepetitions(M_vec)
    outRepet = 0;
    if(length(M_vec)>=5)
        for i=1:(length(M_vec)-2)
            Ai = M_vec(i).A;
            for j=(i+1):(length(M_vec)-1)
                Aj = M_vec(j).A;
                if(sum(Ai==Aj, "all")==64)
                    for k=(j+1):length(M_vec)
                        Ak = M_vec(k).A;
                        if(sum(Ai==Ak, "all")==64)
                            outRepet = 1;
                            return;
                        end
                    end
                    
                end
            end
        end
    end
end

% Evaluate whether it's stalemate for 75 consecutive empty moves
function [outMoves] = isStaleMoves(M_vec)
    outMoves = 0;
    c_u = 0;
    LL = length(M_vec);
    if(LL>150)
        for ii=1:LL
            u=M_vec(ii).captureOrPawnMove;
            if(u==0)
                c_u = c_u+1;
            else
                c_u = 0;
            end
            if(c_u >= 150)
                outMoves = 1;
                return;
            end
        end
    end

end

% Evaluate whether it's check
function [checkStatus] = isCheck(M)
    out_check = 0;
    out_zi = 0;
    out_zj = 0;
    turn = M.turn;
    A = M.A;
    
    sturn = 6+10*turn;
    esci = 0;
    for i=1:8
        for j=1:8
            if(abs(A(i,j)-sturn)<1e-2)
                out_zi = i;
                out_zj = j;
                esci = 1;
                break;
            end
        end
        if(esci==1)
            break;
        end
    end
    if(out_zi>0.5 && out_zj > 0.5)
        for i=1:8
            ic = out_zi + i;
            jc = out_zj;
            if(ic>8.5)
                break;
            end
            if(A(ic,jc) < 0.5)
                continue;
            end
            if(A(ic,jc) == (14-10*turn) || A(ic,jc) == ...
                   (15-10*turn) || (A(ic,jc) == (16-10*turn) && i==1))
                out_check = 1;
                break;
            else
                break;
            end
        end
        if(out_check == 0)
            for i=1:8
                ic = out_zi - i;
                jc = out_zj;
                if(ic<0.5)
                    break;
                end
                if(A(ic,jc) < 0.5)
                    continue;
                end
                if(A(ic,jc) == (14-10*turn) || A(ic,jc) == ...
                       (15-10*turn) || (A(ic,jc) == (16-10*turn) && i==1))
                    out_check = 1;
                    break;
                else
                    break;
                end
            end
        end
        if(out_check == 0)
            for j=1:8
                ic = out_zi;
                jc = out_zj+j;
                if(jc>8.5)
                    break;
                end
                if(A(ic,jc) < 0.5)
                    continue;
                end
                if(A(ic,jc) == (14-10*turn) || A(ic,jc) == ...
                       (15-10*turn) || (A(ic,jc) == (16-10*turn) && j==1))
                    out_check = 1;
                    break;
                else
                    break;
                end
            end
        end
        if(out_check == 0)
            for j=1:8
                ic = out_zi;
                jc = out_zj-j;
                if(jc<0.5)
                    break;
                end
                if(A(ic,jc) < 0.5)
                    continue;
                end
                if(A(ic,jc) == (14-10*turn) || A(ic,jc) == ...
                       (15-10*turn) || (A(ic,jc) == (16-10*turn) && j==1))
                    out_check = 1;
                    break;
                else
                    break;
                end
            end
        end
        if(out_check == 0)
            for k=1:8
                ic = out_zi+k;
                jc = out_zj+k;
                if(ic>8.5 || jc>8.5)
                    break;
                end
                if(A(ic,jc) < 0.5)
                    continue;
                end
                if(A(ic,jc) == (12-10*turn) || A(ic,jc) == ...
                       (15-10*turn) || (A(ic,jc) == (16-10*turn) && k==1))
                    out_check = 1;
                    break;
                else
                    break;
                end
            end
        end
        if(out_check == 0)
            for k=1:8
                ic = out_zi+k;
                jc = out_zj-k;
                if(ic>8.5 || jc<0.5)
                    break;
                end
                if(A(ic,jc) < 0.5)
                    continue;
                end
                if(A(ic,jc) == (12-10*turn) || A(ic,jc) == ...
                       (15-10*turn) || (A(ic,jc) == (16-10*turn) && k==1))
                    out_check = 1;
                    break;
                else
                    break;
                end
            end
        end
        if(out_check == 0)
            for k=1:8
                ic = out_zi-k;
                jc = out_zj+k;
                if(ic<0.5 || jc>8.5)
                    break;
                end
                if(A(ic,jc) < 0.5)
                    continue;
                end
                if(A(ic,jc) == (12-10*turn) || A(ic,jc) == ...
                       (15-10*turn) || (A(ic,jc) == (16-10*turn) && k==1))
                    out_check = 1;
                    break;
                else
                    break;
                end
            end
        end
        if(out_check == 0)
            for k=1:8
                ic = out_zi-k;
                jc = out_zj-k;
                if(ic<0.5 || jc<0.5)
                    break;
                end
                if(A(ic,jc) < 0.5)
                    continue;
                end
                if(A(ic,jc) == (12-10*turn) || A(ic,jc) == ...
                       (15-10*turn) || (A(ic,jc) == (16-10*turn) && k==1))
                    out_check = 1;
                    break;
                else
                    break;
                end
            end
        end
        if(out_check == 0)
            aa = [1, 2,  2,  1, -1, -2, -2, -1];
            bb = [2, 1, -1, -2, -2, -1,  1,  2];
            for k=1:8
                ic = out_zi+aa(k);
                jc = out_zj+bb(k);
                if(ic<0.5 || jc<0.5 || ic>8.5 || jc>8.5)
                    continue;
                end
                if(A(ic,jc) < 0.5)
                    continue;
                end
                if(A(ic,jc) == (13-10*turn))
                    out_check = 1;
                    break;
                else
                    continue;
                end
            end
        end
        if(out_check == 0)
            if(turn==0)
                aa = [-1, -1];
                bb = [1, -1];
            else
                aa = [1, 1];
                bb = [1, -1];
            end
            for k=1:2
                ic = out_zi+aa(k);
                jc = out_zj+bb(k);
                if(ic<0.5 || jc<0.5 || ic>8.5 || jc>8.5)
                    continue;
                end
                if(A(ic,jc) < 0.5)
                    continue;
                end
                if(A(ic,jc) == (11-10*turn))
                    out_check = 1;
                    break;
                else
                    continue;
                end
            end
        end
    end
    checkStatus.out_check = out_check;
    checkStatus.out_zi = out_zi;
    checkStatus.out_zj = out_zj;


end

% Evaluate whether it's checkmate
function [checkmateStatus] = isCheckmate(M)
    checkmateStatus = 0;
    checkStatus = isCheck(M);
    out_check = checkStatus.out_check;
    if(out_check==0)
        return;
    end
    S_out = listNewPossibleConfigurations(M);
    if(isempty(S_out))
        checkmateStatus = 1;
    end

end

% This just calls findBestMove with depths 1, 2, ..., max_depth
function [bestMove, usedDepth] = chooseBestMove(M_vec, max_depth)
    global max_iter
    global max_time
    usedDepth = 0;
    t8=tic();
    for iii=1:max_depth
        max_iter = iii;
        whys=0;
        [bestMove_temp, ~, whys] = findBestMove(M_vec, 0);
        if(whys==1)
            break;
        end
        bestMove = bestMove_temp;
        %bestMove.A
        usedDepth = iii;
        if(toc(t8)>max_time)
            break;
        end
        
    end
    if(iii==1 && whys==1)
        max_iter = 1;
        [bestMove, ~, ~] = findBestMove(M_vec, 0);
        usedDepth = 1;
    end
    max_iter = max_depth;

end

% Longest and most boring function: list of all the available legal moves
function [S_out] = listNewPossibleConfigurations(M, rec0)
    isCheck_res = isCheck(M);
    isCheck_out = isCheck_res.out_check;
    rec1 = 0;
    if(nargin==2)
        rec1 = rec0;
    end
    S_out = [];
    p = 0;
    if(M.turn==0)                                          % WHITE TO MOVE
        for i=1:8
            for j=1:8
                A = M.A;
                Aij = A(i,j);
                if(Aij>0.5 && Aij<10)         
                    if(Aij==1)                            % Pawn movements
                        if(A(i-1,j)==0 && i>2.5)                    % 1 up
                            M_out = M;
                            A_out = A;
                            A_out(i,j) = 0;
                            A_out(i-1,j) = 1;
                            M_out.A = A_out;
                            M_out.lastMove = [i,j,i-1,j];
                            M_out.captureOrPawnMove = 1;
                            salvala = 1;
                            if(rec1==0)
                                checkStatus = isCheck(M_out);
                                if(checkStatus.out_check==1)
                                    salvala = 0;
                                end
                            end
                            if(salvala==1)
                                p = p+1;
                                M_out.turn = 1;
                                S_out(p).M = M_out;
                            end
                        end
                        if(i==7)                                    % 2 up
                            if(A(i-1,j)==0 && A(i-2,j)==0)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i-2,j) = 1;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i-2,j];
                                M_out.captureOrPawnMove = 1;
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                        end
                        if(j>1.5 && i>2.5)                     % Take left
                            if(A(i-1,j-1)>10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i-1,j-1) = 1;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i-1,j-1];
                                M_out.captureOrPawnMove = 1;
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                        end
                        if(j<7.5 && i>2.5)                    % Take right
                            if(A(i-1,j+1)>10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i-1,j+1) = 1;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i-1,j+1];
                                M_out.captureOrPawnMove = 1;
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                        end
                        if(i==2)                               % Promotion
                            if(A(1,j)==0)                             % Up
                                M_out = M;
                                M_out_temp = M_out;
                                M_out_temp.A = A;
                                M_out_temp.A(i,j)=0;
                                M_out_temp.A(1,j)=5;
                                M_out_temp.lastMove = [i,j,1,j];
                                M_out.captureOrPawnMove = 1;
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out_temp);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    M_out.turn = 1;
                                    A_out = A;
                                    A_out(i,j) = 0;
                                    for qq=2:5
                                    A_out(1,j) = qq;
                                    M_out.A = A_out;
                                    M_out.lastMove = [i,j,1,j];
                                    p = p+1;
                                    S_out(p).M = M_out;
                                    end
                                end
                            end
                            if(j>1.5)                          % Take left
                                if(A(1,j-1)>10)
                                    M_out = M;
                                    M_out_temp = M_out;
                                    M_out_temp.A = A;
                                    M_out_temp.A(i,j)=0;
                                    M_out_temp.A(1,j-1)=5;
                                    M_out_temp.lastMove = [i,j,1,j-1];
                                    M_out.captureOrPawnMove = 1;
                                    salvala = 1;
                                    if(rec1==0)
                                        checkStatus = isCheck(M_out_temp);
                                        if(checkStatus.out_check==1)
                                            salvala = 0;
                                        end
                                    end
                                    if(salvala==1)
                                        M_out.turn = 1;
                                        A_out = A;
                                        A_out(i,j) = 0;
                                        for qq=2:5
                                        A_out(1,j-1) = qq;
                                        M_out.A = A_out;
                                        M_out.lastMove = [i,j,1,j-1];
                                        p = p+1;
                                        S_out(p).M = M_out;
                                        end
                                    end
                                end
                            end
                            if(j<7.5)                         % Take right
                                if(A(1,j+1)>10)
                                    M_out = M;
                                    M_out_temp = M_out;
                                    M_out_temp.A = A;
                                    M_out_temp.A(i,j)=0;
                                    M_out_temp.A(1,j+1)=5;
                                    M_out_temp.lastMove = [i,j,1,j+1];
                                    M_out.captureOrPawnMove = 1;
                                    salvala = 1;
                                    if(rec1==0)
                                        checkStatus = isCheck(M_out_temp);
                                        if(checkStatus.out_check==1)
                                            salvala = 0;
                                        end
                                    end
                                    if(salvala==1)
                                        M_out.turn = 1;
                                        A_out = A;
                                        A_out(i,j) = 0;
                                        for qq=2:5
                                        A_out(1,j+1) = qq;
                                        M_out.A = A_out;
                                        M_out.lastMove = [i,j,1,j+1];
                                        p = p+1;
                                        S_out(p).M = M_out;
                                        end
                                    end
                                end
                            end
                        end
                        if(i==4)                              % En passant
                            if(j>1.5)                               % Left
                                if(A(3,j-1)==0 && A(4,j-1)==11 &&...
                                        (sum(M.lastMove ==...
                                        [2,j-1,4,j-1])==4))
                                    M_out = M;
                                    A_out = A;
                                    A_out(i,j) = 0;
                                    A_out(4,j-1) = 0;
                                    A_out(3,j-1) = 1;
                                    M_out.A = A_out;
                                    M_out.lastMove = [i,j,3,j-1];
                                    M_out.captureOrPawnMove = 1;
                                    salvala = 1;
                                    if(rec1==0)
                                        checkStatus = isCheck(M_out);
                                        if(checkStatus.out_check==1)
                                            salvala = 0;
                                        end
                                    end
                                    if(salvala==1)
                                        p = p+1;
                                        M_out.turn = 1;
                                        S_out(p).M = M_out;
                                    end
                                end

                            end
                            if(j<7.5)                              % Right
                                if(A(3,j+1)==0 && A(4,j+1)==11 &&...
                                        (sum(M.lastMove ==...
                                        [2,j+1,4,j+1])==4))
                                    M_out = M;
                                    A_out = A;
                                    A_out(i,j) = 0;
                                    A_out(4,j+1) = 0;
                                    A_out(3,j+1) = 1;
                                    M_out.A = A_out;
                                    M_out.lastMove = [i,j,3,j+1];
                                    M_out.captureOrPawnMove = 1;
                                    salvala = 1;
                                    if(rec1==0)
                                        checkStatus = isCheck(M_out);
                                        if(checkStatus.out_check==1)
                                            salvala = 0;
                                        end
                                    end
                                    if(salvala==1)
                                        p = p+1;
                                        M_out.turn = 1;
                                        S_out(p).M = M_out;
                                    end
                                end

                            end

                        end
    
                    elseif(Aij==2)                      % Bishop movements
                        for ijk=(1:1:7)
                            notAdm = ((i+ijk)>8.5) || ((j+ijk)>8.5);
                            if(notAdm==1)
                                break;
                            end
                            if(A(i+ijk,j+ijk)==0||A(i+ijk,j+ijk)>10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+ijk,j+ijk) = 2;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+ijk,j+ijk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i+ijk,j+ijk)>10)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i+ijk,j+ijk)>0.5)
                                break;
                            end
                        end
                        for ijk=(1:1:7)
                            notAdm = ((i+ijk)>8.5) || ((j-ijk)<0.5);
                            if(notAdm==1)
                                break;
                            end
                            if(A(i+ijk,j-ijk)==0||A(i+ijk,j-ijk)>10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+ijk,j-ijk) = 2;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+ijk,j-ijk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i+ijk,j-ijk)>10)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i+ijk,j-ijk)>0.5)
                                break;
                            end
                        end
                        for ijk=(1:1:7)
                            notAdm = ((i-ijk)<0.5) || ((j+ijk)>8.5);
                            if(notAdm==1)
                                break;
                            end
                            if(A(i-ijk,j+ijk)==0||A(i-ijk,j+ijk)>10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i-ijk,j+ijk) = 2;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i-ijk,j+ijk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i-ijk,j+ijk)>10)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i-ijk,j+ijk)>0.5)
                                break;
                            end
                        end
                        for ijk=(1:1:7)
                            notAdm = ((i-ijk)<0.5) || ((j-ijk)<0.5);
                            if(notAdm==1)
                                break;
                            end
                            if(A(i-ijk,j-ijk)==0||A(i-ijk,j-ijk)>10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i-ijk,j-ijk) = 2;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i-ijk,j-ijk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i-ijk,j-ijk)>10)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i-ijk,j-ijk)>0.5)
                                break;
                            end
                        end

                    elseif(Aij==3)                      % kNight movements
                        for ik=1:8
                            for jk=1:8
                                if((A(ik,jk)==0 || A(ik,jk)>10) &&...
                                        abs(norm([ik-i,jk-j])-...
                                        sqrt(5))<1e-2)
                                    M_out = M;
                                    A_out = A;
                                    A_out(i,j) = 0;
                                    A_out(ik,jk) = 3;
                                    M_out.A = A_out;
                                    M_out.lastMove = [i,j,ik,jk];
                                    M_out.captureOrPawnMove = 0;
                                    if(A(ik,jk)>10)
                                        M_out.captureOrPawnMove = 1;
                                    end
                                    salvala = 1;
                                    if(rec1==0)
                                        checkStatus = isCheck(M_out);
                                        if(checkStatus.out_check==1)
                                            salvala = 0;
                                        end
                                    end
                                    if(salvala==1)
                                        p = p+1;
                                        M_out.turn = 1;
                                        S_out(p).M = M_out;
                                    end
                                end

                            end
                        end
                    elseif(Aij==4)                        % Rook movements
                        for ik=(1:1:7)                                % Up
                            notAdm = (i+ik)>8.5;
                            if(notAdm==1)
                                break;
                            end
                            if(A(i+ik,j)==0||A(i+ik,j)>10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+ik,j) = 4;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+ik,j];
                                M_out.captureOrPawnMove = 0;
                                if(A(i+ik,j)>10)
                                    M_out.captureOrPawnMove = 1;
                                end
                                if(i==8&&j==1)
                                    M_out.whiteLeftRookMoved = 1;
                                elseif(i==8&&j==8)
                                    M_out.whiteRightRookMoved = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i+ik,j)>0.5)
                                break;
                            end
                        end
                        for ik=(-1:-1:-7)                           % Down
                            notAdm = (i+ik)<0.5;
                            if(notAdm==1)
                                break;
                            end
                            if(A(i+ik,j)==0||A(i+ik,j)>10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+ik,j) = 4;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+ik,j];
                                M_out.captureOrPawnMove = 0;
                                if(A(i+ik,j)>10)
                                    M_out.captureOrPawnMove = 1;
                                end
                                if(i==8&&j==1)
                                    M_out.whiteLeftRookMoved = 1;
                                elseif(i==8&&j==8)
                                    M_out.whiteRightRookMoved = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i+ik,j)>0.5)
                                break;
                            end
                        end
                        for jk=(1:1:7)                             % Right
                            notAdm = (j+jk)>8.5;
                            if(notAdm==1)
                                break;
                            end
                            if(A(i,j+jk)==0||A(i,j+jk)>10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i,j+jk) = 4;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i,j+jk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i,j+jk)>10)
                                    M_out.captureOrPawnMove = 1;
                                end
                                if(i==8&&j==1)
                                    M_out.whiteLeftRookMoved = 1;
                                elseif(i==8&&j==8)
                                    M_out.whiteRightRookMoved = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i,j+jk)>0.5)
                                break;
                            end
                        end
                        for jk=(-1:-1:-7)                           % Left
                            notAdm = (j+jk)<0.5;
                            if(notAdm==1)
                                break;
                            end
                            if(A(i,j+jk)==0||A(i,j+jk)>10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i,j+jk) = 4;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i,j+jk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i,j+jk)>10)
                                    M_out.captureOrPawnMove = 1;
                                end
                                if(i==8&&j==1)
                                    M_out.whiteLeftRookMoved = 1;
                                elseif(i==8&&j==8)
                                    M_out.whiteRightRookMoved = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i,j+jk)>0.5)
                                break;
                            end
                        end

                    elseif(Aij==5)                       % Queen movements
                        for ijk=(1:1:7)
                            notAdm = ((i+ijk)>8.5) || ((j+ijk)>8.5);
                            if(notAdm==1)
                                break;
                            end
                            if(A(i+ijk,j+ijk)==0||A(i+ijk,j+ijk)>10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+ijk,j+ijk) = 5;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+ijk,j+ijk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i+ijk,j+ijk)>10)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i+ijk,j+ijk)>0.5)
                                break;
                            end
                        end
                        for ijk=(1:1:7)
                            notAdm = ((i+ijk)>8.5) || ((j-ijk)<0.5);
                            if(notAdm==1)
                                break;
                            end
                            if(A(i+ijk,j-ijk)==0||A(i+ijk,j-ijk)>10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+ijk,j-ijk) = 5;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+ijk,j-ijk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i+ijk,j-ijk)>10)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i+ijk,j-ijk)>0.5)
                                break;
                            end
                        end
                        for ijk=(1:1:7)
                            notAdm = ((i-ijk)<0.5) || ((j+ijk)>8.5);
                            if(notAdm==1)
                                break;
                            end
                            if(A(i-ijk,j+ijk)==0||A(i-ijk,j+ijk)>10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i-ijk,j+ijk) = 5;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i-ijk,j+ijk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i-ijk,j+ijk)>10)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i-ijk,j+ijk)>0.5)
                                break;
                            end
                        end
                        for ijk=(1:1:7)
                            notAdm = ((i-ijk)<0.5) || ((j-ijk)<0.5);
                            if(notAdm==1)
                                break;
                            end
                            if(A(i-ijk,j-ijk)==0||A(i-ijk,j-ijk)>10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i-ijk,j-ijk) = 5;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i-ijk,j-ijk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i-ijk,j-ijk)>10)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i-ijk,j-ijk)>0.5)
                                break;
                            end
                        end
                        for ik=(1:1:7)                                % Up
                            notAdm = (i+ik)>8.5;
                            if(notAdm==1)
                                break;
                            end
                            if(A(i+ik,j)==0||A(i+ik,j)>10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+ik,j) = 5;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+ik,j];
                                M_out.captureOrPawnMove = 0;
                                if(A(i+ik,j)>10)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i+ik,j)>0.5)
                                break;
                            end
                        end
                        for ik=(-1:-1:-7)                           % Down
                            notAdm = (i+ik)<0.5;
                            if(notAdm==1)
                                break;
                            end
                            if(A(i+ik,j)==0||A(i+ik,j)>10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+ik,j) = 5;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+ik,j];
                                M_out.captureOrPawnMove = 0;
                                if(A(i+ik,j)>10)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i+ik,j)>0.5)
                                break;
                            end
                        end
                        for jk=(1:1:7)                             % Right
                            notAdm = (j+jk)>8.5;
                            if(notAdm==1)
                                break;
                            end
                            if(A(i,j+jk)==0||A(i,j+jk)>10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i,j+jk) = 5;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i,j+jk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i,j+jk)>10)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i,j+jk)>0.5)
                                break;
                            end
                        end
                        for jk=(-1:-1:-7)                           % Left
                            notAdm = (j+jk)<0.5;
                            if(notAdm==1)
                                break;
                            end
                            if(A(i,j+jk)==0||A(i,j+jk)>10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i,j+jk) = 5;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i,j+jk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i,j+jk)>10)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 1;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i,j+jk)>0.5)
                                break;
                            end
                        end

                    elseif(Aij==6)                        % King movements
                        for ik=-1:1:1
                            for jk=-1:1:1
                                notAdm = ((i+ik)<0.5)||((i+ik)>8.5)||...
                                    ((j+jk)<0.5)||((j+jk)>8.5);
                                if((ik==0 && jk==0)||notAdm==1)
                                    continue;
                                end
                                if(A(i+ik,j+jk)==0||A(i+ik,j+jk)>10)
                                    M_out = M;
                                    A_out = A;
                                    A_out(i,j) = 0;
                                    A_out(i+ik,j+jk) = 6;
                                    M_out.A = A_out;
                                    M_out.lastMove = [i,j,i+ik,j+jk];
                                    M_out.captureOrPawnMove = 0;
                                    if(A(i+ik,j+jk)>10)
                                        M_out.captureOrPawnMove = 1;
                                    end
                                    M_out.whiteKingMoved = 1;
                                    salvala = 1;
                                    if(rec1==0)
                                        checkStatus = isCheck(M_out);
                                        if(checkStatus.out_check==1)
                                            salvala = 0;
                                        end
                                    end
                                    if(salvala==1)
                                        p = p+1;
                                        M_out.turn = 1;
                                        S_out(p).M = M_out;
                                    end
                                end
                            end
                        end
                        if(M.whiteAlreadyCastled==0 && ...
                                M.whiteKingMoved==0 && ...
                                isCheck_out == 0)                 % Castle
                            if(M.whiteLeftRookMoved==0)             % Long
                                if(A(8,1)==4&&A(8,2)==0&&...
                                        A(8,3)==0&&A(8,4)==0)
                                    M_out = M;
                                    A_out = A;
                                    A_out(i,j) = 0;
                                    A_out(8,1) = 0;
                                    A_out(8,2) = 0;
                                    A_out(8,3) = 6;
                                    A_out(8,4) = 4;
                                    M_out.A = A_out;
                                    M_out.lastMove = [i,j,8,3];
                                    M_out.captureOrPawnMove = 0;
                                    M_out.whiteKingMoved = 1;
                                    M_out.whiteAlreadyCastled = 1;
                                    salvala = 1;
                                    if(rec1==0)
                                        M_out_temp = M_out;
                                        M_out_temp.A(8,1) = 4;
                                        M_out_temp.A(8,2) = 0;
                                        M_out_temp.A(8,3) = 0;
                                        M_out_temp.A(8,4) = 6;
                                        M_out_temp.A(8,5) = 0;
                                        checkStatus = isCheck(M_out_temp);
                                        if(checkStatus.out_check==1)
                                            salvala = 0;
                                        end
                                        if(salvala == 1)
                                            M_out_temp = M_out;
                                            M_out_temp.A(8,1) = 4;
                                            M_out_temp.A(8,2) = 0;
                                            M_out_temp.A(8,3) = 6;
                                            M_out_temp.A(8,4) = 0;
                                            M_out_temp.A(8,5) = 0;
                                            checkStatus =...
                                                isCheck(M_out_temp);
                                            if(checkStatus.out_check==1)
                                                salvala = 0;
                                            end
                                        end
                                        if(salvala == 1)
                                            M_out_temp = M_out;
                                            checkStatus =...
                                                isCheck(M_out_temp);
                                            if(checkStatus.out_check==1)
                                                salvala = 0;
                                            end
                                        end
                                    end
                                    if(salvala==1)
                                        p = p+1;
                                        M_out.turn = 1;
                                        S_out(p).M = M_out;
                                    end
                                end
                            end
                            if(M.whiteRightRookMoved==0)           % Short
                                if(A(8,8)==4&&A(8,7)==0&&A(8,6)==0)
                                    M_out = M;
                                    A_out = A;
                                    A_out(i,j) = 0;
                                    A_out(8,8) = 0;
                                    A_out(8,7) = 6;
                                    A_out(8,6) = 4;
                                    M_out.A = A_out;
                                    M_out.lastMove = [i,j,8,7];
                                    M_out.captureOrPawnMove = 0;
                                    M_out.whiteKingMoved = 1;
                                    M_out.whiteAlreadyCastled = 1;
                                    salvala = 1;
                                    if(rec1==0)
                                        M_out_temp = M_out;
                                        M_out_temp.A(8,5) = 0;
                                        M_out_temp.A(8,6) = 6;
                                        M_out_temp.A(8,7) = 0;
                                        M_out_temp.A(8,8) = 4;
                                        checkStatus = isCheck(M_out_temp);
                                        if(checkStatus.out_check==1)
                                            salvala = 0;
                                        end
                                        if(salvala == 1)
                                            M_out_temp = M_out;
                                            M_out_temp.A(8,5) = 0;
                                            M_out_temp.A(8,6) = 0;
                                            M_out_temp.A(8,7) = 6;
                                            M_out_temp.A(8,8) = 4;
                                            checkStatus =...
                                                isCheck(M_out_temp);
                                            if(checkStatus.out_check==1)
                                                salvala = 0;
                                            end
                                        end
                                        if(salvala == 1)
                                            M_out_temp = M_out;
                                            checkStatus =...
                                                isCheck(M_out_temp);
                                            if(checkStatus.out_check==1)
                                                salvala = 0;
                                            end
                                        end
                                    end
                                    if(salvala==1)
                                        p = p+1;
                                        M_out.turn = 1;
                                        S_out(p).M = M_out;
                                    end
                                end
                            end

                        end

                    end
                end
            end
    
        end
    elseif(M.turn==1)                                      % BLACK TO MOVE
        for i=1:8
            for j=1:8
                A = M.A;
                Aij = A(i,j);
                if(Aij>10)         
                    if(Aij==11)                           % Pawn movements
                        if(A(i+1,j)==0 && i<6.5)                    % 1 up
                            M_out = M;
                            A_out = A;
                            A_out(i,j) = 0;
                            A_out(i+1,j) = 11;
                            M_out.A = A_out;
                            M_out.lastMove = [i,j,i+1,j];
                            M_out.captureOrPawnMove = 1;
                            salvala = 1;
                            if(rec1==0)
                                checkStatus = isCheck(M_out);
                                if(checkStatus.out_check==1)
                                    salvala = 0;
                                end
                            end
                            if(salvala==1)
                                p = p+1;
                                M_out.turn = 0;
                                S_out(p).M = M_out;
                            end
                        end
                        if(i==2)                                    % 2 up
                            if(A(i+1,j)==0 && A(i+2,j)==0)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+2,j) = 11;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+2,j];
                                M_out.captureOrPawnMove = 1;
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                        end
                        if(j>1.5 && i<6.5)                     % Take left
                            if(A(i+1,j-1)>0.5 && A(i+1,j-1)<10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+1,j-1) = 11;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+1,j-1];
                                M_out.captureOrPawnMove = 1;
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                        end
                        if(j<7.5 && i<6.5)                    % Take right
                            if(A(i+1,j+1)>0.5 && A(i+1,j+1)<10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+1,j+1) = 11;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+1,j+1];
                                M_out.captureOrPawnMove = 1;
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                        end
                        if(i==7)                               % Promotion
                            if(A(8,j)==0)                             % Up
                                M_out = M;
                                M_out_temp = M_out;
                                M_out_temp.A = A;
                                M_out_temp.A(i,j)=0;
                                M_out_temp.A(8,j)=15;
                                M_out_temp.lastMove = [i,j,8,j];
                                M_out.captureOrPawnMove = 1;
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out_temp);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    M_out.turn = 0;
                                    A_out = A;
                                    A_out(i,j) = 0;
                                    for qq=12:15
                                    A_out(8,j) = qq;
                                    M_out.A = A_out;
                                    M_out.lastMove = [i,j,8,j];
                                    p = p+1;
                                    S_out(p).M = M_out;
                                    end
                                end
                            end
                            if(j>1.5)                          % Take left
                                if(A(8,j-1)>0.5 && A(8,j-1)<10)
                                    M_out = M;
                                    M_out_temp = M_out;
                                    M_out_temp.A = A;
                                    M_out_temp.A(i,j)=0;
                                    M_out_temp.A(8,j-1)=15;
                                    M_out_temp.lastMove = [i,j,8,j-1];
                                    M_out.captureOrPawnMove = 1;
                                    salvala = 1;
                                    if(rec1==0)
                                        checkStatus = isCheck(M_out_temp);
                                        if(checkStatus.out_check==1)
                                            salvala = 0;
                                        end
                                    end
                                    if(salvala==1)
                                        M_out.turn = 0;
                                        A_out = A;
                                        A_out(i,j) = 0;
                                        for qq=12:15
                                        A_out(8,j-1) = qq;
                                        M_out.A = A_out;
                                        M_out.lastMove = [i,j,8,j-1];
                                        p = p+1;
                                        S_out(p).M = M_out;
                                        end
                                    end
                                end
                            end
                            if(j<7.5)                         % Take right
                                if(A(8,j+1)>0.5 && A(8,j+1)<10)
                                    M_out = M;
                                    M_out_temp = M_out;
                                    M_out_temp.A = A;
                                    M_out_temp.A(i,j)=0;
                                    M_out_temp.A(8,j+1)=15;
                                    M_out_temp.lastMove = [i,j,8,j+1];
                                    M_out.captureOrPawnMove = 1;
                                    salvala = 1;
                                    if(rec1==0)
                                        checkStatus = isCheck(M_out_temp);
                                        if(checkStatus.out_check==1)
                                            salvala = 0;
                                        end
                                    end
                                    if(salvala==1)
                                        M_out.turn = 0;
                                        A_out = A;
                                        A_out(i,j) = 0;
                                        for qq=12:15
                                        A_out(8,j+1) = qq;
                                        M_out.A = A_out;
                                        M_out.lastMove = [i,j,8,j+1];
                                        p = p+1;
                                        S_out(p).M = M_out;
                                        end
                                    end
                                end
                            end
                        end
                        if(i==5)                              % En passant
                            if(j>1.5)                               % Left
                                if(A(6,j-1)==0 && A(5,j-1)==1 &&...
                                        (sum(M.lastMove ==...
                                        [7,j-1,5,j-1])==4))
                                    M_out = M;
                                    A_out = A;
                                    A_out(i,j) = 0;
                                    A_out(5,j-1) = 0;
                                    A_out(6,j-1) = 11;
                                    M_out.A = A_out;
                                    M_out.lastMove = [i,j,6,j-1];
                                    M_out.captureOrPawnMove = 1;
                                    salvala = 1;
                                    if(rec1==0)
                                        checkStatus = isCheck(M_out);
                                        if(checkStatus.out_check==1)
                                            salvala = 0;
                                        end
                                    end
                                    if(salvala==1)
                                        p = p+1;
                                        M_out.turn = 0;
                                        S_out(p).M = M_out;
                                    end
                                end

                            end
                            if(j<7.5)                              % Right
                                if(A(6,j+1)==0 && A(5,j+1)==11 && ...
                                        (sum(M.lastMove ==...
                                        [7,j+1,5,j+1])==4))
                                    M_out = M;
                                    A_out = A;
                                    A_out(i,j) = 0;
                                    A_out(5,j+1) = 0;
                                    A_out(6,j+1) = 11;
                                    M_out.A = A_out;
                                    M_out.lastMove = [i,j,6,j+1];
                                    M_out.captureOrPawnMove = 1;
                                    salvala = 1;
                                    if(rec1==0)
                                        checkStatus = isCheck(M_out);
                                        if(checkStatus.out_check==1)
                                            salvala = 0;
                                        end
                                    end
                                    if(salvala==1)
                                        p = p+1;
                                        M_out.turn = 0;
                                        S_out(p).M = M_out;
                                    end
                                end

                            end

                        end
    
                    elseif(Aij==12)                     % Bishop movements
                        for ijk=(1:1:7)
                            notAdm = ((i+ijk)>8.5) || ((j+ijk)>8.5);
                            if(notAdm==1)
                                break;
                            end
                            if(A(i+ijk,j+ijk)==0||A(i+ijk,j+ijk)<10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+ijk,j+ijk) = 12;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+ijk,j+ijk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i+ijk,j+ijk)>1)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i+ijk,j+ijk)>0.5)
                                break;
                            end
                        end
                        for ijk=(1:1:7)
                            notAdm = ((i+ijk)>8.5) || ((j-ijk)<0.5);
                            if(notAdm==1)
                                break;
                            end
                            if(A(i+ijk,j-ijk)==0||A(i+ijk,j-ijk)<10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+ijk,j-ijk) = 12;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+ijk,j-ijk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i+ijk,j-ijk)>1)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i+ijk,j-ijk)>0.5)
                                break;
                            end
                        end
                        for ijk=(1:1:7)
                            notAdm = ((i-ijk)<0.5) || ((j+ijk)>8.5);
                            if(notAdm==1)
                                break;
                            end
                            if(A(i-ijk,j+ijk)==0||A(i-ijk,j+ijk)<10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i-ijk,j+ijk) = 12;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i-ijk,j+ijk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i-ijk,j+ijk)>1)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i-ijk,j+ijk)>0.5)
                                break;
                            end
                        end
                        for ijk=(1:1:7)
                            notAdm = ((i-ijk)<0.5) || ((j-ijk)<0.5);
                            if(notAdm==1)
                                break;
                            end
                            if(A(i-ijk,j-ijk)==0||A(i-ijk,j-ijk)<10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i-ijk,j-ijk) = 12;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i-ijk,j-ijk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i-ijk,j-ijk)>1)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i-ijk,j-ijk)>0.5)
                                break;
                            end
                        end

                    elseif(Aij==13)                     % kNight movements
                        for ik=1:8
                            for jk=1:8
                                if((A(ik,jk)==0 || A(ik,jk)<10) &&...
                                      abs(norm([ik-i,jk-j])-sqrt(5))<1e-2)
                                    M_out = M;
                                    A_out = A;
                                    A_out(i,j) = 0;
                                    A_out(ik,jk) = 13;
                                    M_out.A = A_out;
                                    M_out.lastMove = [i,j,ik,jk];
                                    M_out.captureOrPawnMove = 0;
                                    if(A(ik,jk)>1)
                                        M_out.captureOrPawnMove = 1;
                                    end
                                    salvala = 1;
                                    if(rec1==0)
                                        checkStatus = isCheck(M_out);
                                        if(checkStatus.out_check==1)
                                            salvala = 0;
                                        end
                                    end
                                    if(salvala==1)
                                        p = p+1;
                                        M_out.turn = 0;
                                        S_out(p).M = M_out;
                                    end
                                end

                            end
                        end
                    elseif(Aij==14)                       % Rook movements
                        for ik=(1:1:7)                                % Up
                            notAdm = (i+ik)>8.5;
                            if(notAdm==1)
                                break;
                            end
                            if(A(i+ik,j)==0||A(i+ik,j)<10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+ik,j) = 14;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+ik,j];
                                M_out.captureOrPawnMove = 0;
                                if(A(i+ik,j)>1)
                                    M_out.captureOrPawnMove = 1;
                                end
                                if(i==1&&j==1)
                                    M_out.blackLeftRookMoved = 1;
                                elseif(i==1&&j==8)
                                    M_out.blackRightRookMoved = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i+ik,j)>0.5)
                                break;
                            end
                        end
                        for ik=(-1:-1:-7)                           % Down
                            notAdm = (i+ik)<0.5;
                            if(notAdm==1)
                                break;
                            end
                            if(A(i+ik,j)==0||A(i+ik,j)<10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+ik,j) = 14;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+ik,j];
                                M_out.captureOrPawnMove = 0;
                                if(A(i+ik,j)>1)
                                    M_out.captureOrPawnMove = 1;
                                end
                                if(i==1&&j==1)
                                    M_out.blackLeftRookMoved = 1;
                                elseif(i==1&&j==8)
                                    M_out.blackRightRookMoved = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i+ik,j)>0.5)
                                break;
                            end
                        end
                        for jk=(1:1:7)                             % Right
                            notAdm = (j+jk)>8.5;
                            if(notAdm==1)
                                break;
                            end
                            if(A(i,j+jk)==0||A(i,j+jk)<10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i,j+jk) = 14;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i,j+jk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i,j+jk)>1)
                                    M_out.captureOrPawnMove = 1;
                                end
                                if(i==1&&j==1)
                                    M_out.blackLeftRookMoved = 1;
                                elseif(i==1&&j==8)
                                    M_out.blackRightRookMoved = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i,j+jk)>0.5)
                                break;
                            end
                        end
                        for jk=(-1:-1:-7)                           % Left
                            notAdm = (j+jk)<0.5;
                            if(notAdm==1)
                                break;
                            end
                            if(A(i,j+jk)==0||A(i,j+jk)<10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i,j+jk) = 14;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i,j+jk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i,j+jk)>1)
                                    M_out.captureOrPawnMove = 1;
                                end
                                if(i==1&&j==1)
                                    M_out.blackLeftRookMoved = 1;
                                elseif(i==1&&j==8)
                                    M_out.blackRightRookMoved = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i,j+jk)>0.5)
                                break;
                            end
                        end

                    elseif(Aij==15)                      % Queen movements
                        for ijk=(1:1:7)
                            notAdm = ((i+ijk)>8.5) || ((j+ijk)>8.5);
                            if(notAdm==1)
                                break;
                            end
                            if(A(i+ijk,j+ijk)==0||A(i+ijk,j+ijk)<10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+ijk,j+ijk) = 15;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+ijk,j+ijk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i+ijk,j+ijk)>1)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i+ijk,j+ijk)>0.5)
                                break;
                            end
                        end
                        for ijk=(1:1:7)
                            notAdm = ((i+ijk)>8.5) || ((j-ijk)<0.5);
                            if(notAdm==1)
                                break;
                            end
                            if(A(i+ijk,j-ijk)==0||A(i+ijk,j-ijk)<10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+ijk,j-ijk) = 15;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+ijk,j-ijk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i+ijk,j-ijk)>1)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i+ijk,j-ijk)>0.5)
                                break;
                            end
                        end
                        for ijk=(1:1:7)
                            notAdm = ((i-ijk)<0.5) || ((j+ijk)>8.5);
                            if(notAdm==1)
                                break;
                            end
                            if(A(i-ijk,j+ijk)==0||A(i-ijk,j+ijk)<10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i-ijk,j+ijk) = 15;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i-ijk,j+ijk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i-ijk,j+ijk)>1)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i-ijk,j+ijk)>0.5)
                                break;
                            end
                        end
                        for ijk=(1:1:7)
                            notAdm = ((i-ijk)<0.5) || ((j-ijk)<0.5);
                            if(notAdm==1)
                                break;
                            end
                            if(A(i-ijk,j-ijk)==0||A(i-ijk,j-ijk)<10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i-ijk,j-ijk) = 15;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i-ijk,j-ijk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i-ijk,j-ijk)>1)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i-ijk,j-ijk)>0.5)
                                break;
                            end
                        end
                        for ik=(1:1:7)                                % Up
                            notAdm = (i+ik)>8.5;
                            if(notAdm==1)
                                break;
                            end
                            if(A(i+ik,j)==0||A(i+ik,j)<10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+ik,j) = 15;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+ik,j];
                                M_out.captureOrPawnMove = 0;
                                if(A(i+ik,j)>1)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i+ik,j)>0.5)
                                break;
                            end
                        end
                        for ik=(-1:-1:-7)                           % Down
                            notAdm = (i+ik)<0.5;
                            if(notAdm==1)
                                break;
                            end
                            if(A(i+ik,j)==0||A(i+ik,j)<10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i+ik,j) = 15;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i+ik,j];
                                M_out.captureOrPawnMove = 0;
                                if(A(i+ik,j)>1)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i+ik,j)>0.5)
                                break;
                            end
                        end
                        for jk=(1:1:7)                             % Right
                            notAdm = (j+jk)>8.5;
                            if(notAdm==1)
                                break;
                            end
                            if(A(i,j+jk)==0||A(i,j+jk)<10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i,j+jk) = 15;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i,j+jk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i,j+jk)>1)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i,j+jk)>0.5)
                                break;
                            end
                        end
                        for jk=(-1:-1:-7)                           % Left
                            notAdm = (j+jk)<0.5;
                            if(notAdm==1)
                                break;
                            end
                            if(A(i,j+jk)==0||A(i,j+jk)<10)
                                M_out = M;
                                A_out = A;
                                A_out(i,j) = 0;
                                A_out(i,j+jk) = 15;
                                M_out.A = A_out;
                                M_out.lastMove = [i,j,i,j+jk];
                                M_out.captureOrPawnMove = 0;
                                if(A(i,j+jk)>1)
                                    M_out.captureOrPawnMove = 1;
                                end
                                salvala = 1;
                                if(rec1==0)
                                    checkStatus = isCheck(M_out);
                                    if(checkStatus.out_check==1)
                                        salvala = 0;
                                    end
                                end
                                if(salvala==1)
                                    p = p+1;
                                    M_out.turn = 0;
                                    S_out(p).M = M_out;
                                end
                            end
                            if(A(i,j+jk)>0.5)
                                break;
                            end
                        end

                    elseif(Aij==16)                       % King movements
                        for ik=-1:1:1
                            for jk=-1:1:1
                                notAdm = ((i+ik)<0.5)||((i+ik)>8.5)||...
                                    ((j+jk)<0.5)||((j+jk)>8.5);
                                if((ik==0 && jk==0)||notAdm==1)
                                    continue;
                                end
                                if(A(i+ik,j+jk)==0||A(i+ik,j+jk)<10)
                                    M_out = M;
                                    A_out = A;
                                    A_out(i,j) = 0;
                                    A_out(i+ik,j+jk) = 16;
                                    M_out.A = A_out;
                                    M_out.lastMove = [i,j,i+ik,j+jk];
                                    M_out.captureOrPawnMove = 0;
                                    if(A(i+ik,j+jk)>1)
                                        M_out.captureOrPawnMove = 1;
                                    end
                                    M_out.blackKingMoved = 1;
                                    salvala = 1;
                                    if(rec1==0)
                                        checkStatus = isCheck(M_out);
                                        if(checkStatus.out_check==1)
                                            salvala = 0;
                                        end
                                    end
                                    if(salvala==1)
                                        p = p+1;
                                        M_out.turn = 0;
                                        S_out(p).M = M_out;
                                    end
                                end
                            end
                        end
                        if(M.blackAlreadyCastled==0 &&...
                                M.blackKingMoved==0 && ...
                                isCheck_out==0)                   % Castle
                            if(M.blackLeftRookMoved==0)      % Long (left)
                                if(A(1,1)==14&&A(1,2)==0&&...
                                        A(1,3)==0&&A(1,4)==0)
                                    M_out = M;
                                    A_out = A;
                                    A_out(i,j) = 0;
                                    A_out(1,1) = 0;
                                    A_out(1,2) = 0;
                                    A_out(1,3) = 16;
                                    A_out(1,4) = 14;
                                    M_out.A = A_out;
                                    M_out.lastMove = [i,j,1,3];
                                    M_out.captureOrPawnMove = 0;
                                    M_out.blackKingMoved = 1;
                                    M_out.blackAlreadyCastled = 1;
                                    salvala = 1;
                                    if(rec1==0)
                                        M_out_temp = M_out;
                                        M_out_temp.A(1,1) = 4;
                                        M_out_temp.A(1,2) = 0;
                                        M_out_temp.A(1,3) = 0;
                                        M_out_temp.A(1,4) = 6;
                                        M_out_temp.A(1,5) = 0;
                                        checkStatus = isCheck(M_out_temp);
                                        if(checkStatus.out_check==1)
                                            salvala = 0;
                                        end
                                        if(salvala == 1)
                                            M_out_temp = M_out;
                                            M_out_temp.A(1,1) = 4;
                                            M_out_temp.A(1,2) = 0;
                                            M_out_temp.A(1,3) = 6;
                                            M_out_temp.A(1,4) = 0;
                                            M_out_temp.A(1,5) = 0;
                                            checkStatus =...
                                                isCheck(M_out_temp);
                                            if(checkStatus.out_check==1)
                                                salvala = 0;
                                            end
                                        end
                                        if(salvala == 1)
                                            M_out_temp = M_out;
                                            checkStatus =...
                                                isCheck(M_out_temp);
                                            if(checkStatus.out_check==1)
                                                salvala = 0;
                                            end
                                        end
                                    end
                                    if(salvala==1)
                                        p = p+1;
                                        M_out.turn = 0;
                                        S_out(p).M = M_out;
                                    end
                                end
                            end
                            if(M.blackRightRookMoved==0)   % Short (right)
                                if(A(1,8)==14&&A(1,7)==0&&A(1,6)==0)
                                    M_out = M;
                                    A_out = A;
                                    A_out(i,j) = 0;
                                    A_out(1,8) = 0;
                                    A_out(1,7) = 16;
                                    A_out(1,6) = 14;
                                    M_out.A = A_out;
                                    M_out.lastMove = [i,j,1,7];
                                    M_out.captureOrPawnMove = 0;
                                    M_out.blackKingMoved = 1;
                                    M_out.blackAlreadyCastled = 1;
                                    salvala = 1;
                                    if(rec1==0)
                                        M_out_temp = M_out;
                                        M_out_temp.A(1,5) = 0;
                                        M_out_temp.A(1,6) = 6;
                                        M_out_temp.A(1,7) = 0;
                                        M_out_temp.A(1,8) = 4;
                                        checkStatus = isCheck(M_out_temp);
                                        if(checkStatus.out_check==1)
                                            salvala = 0;
                                        end
                                        if(salvala == 1)
                                            M_out_temp = M_out;
                                            M_out_temp.A(1,5) = 0;
                                            M_out_temp.A(1,6) = 0;
                                            M_out_temp.A(1,7) = 6;
                                            M_out_temp.A(1,8) = 4;
                                            checkStatus =...
                                                isCheck(M_out_temp);
                                            if(checkStatus.out_check==1)
                                                salvala = 0;
                                            end
                                        end
                                        if(salvala == 1)
                                            M_out_temp = M_out;
                                            checkStatus =...
                                                isCheck(M_out_temp);
                                            if(checkStatus.out_check==1)
                                                salvala = 0;
                                            end
                                        end
                                    end
                                    if(salvala==1)
                                        p = p+1;
                                        M_out.turn = 0;
                                        S_out(p).M = M_out;
                                    end
                                end
                            end

                        end

                    end
                end
            end
    
        end
    end

end

%% Interface
% Plot the base board
function [] = plotConfig(M)
    A = M.A;
    checkStatus = M.checkStatus;
    checkmateStatus = M.checkmateStatus;
    staleStatus = M.staleStatus;
    global whiteBackgroundColor
    global blackBackgroundColor
    global playAlone
    global mvv
    global prev
    global sside
    cla
    f1 = figure(1);
    hold on
    f1.Units = 'pixels';
    f1.Position = 0.8*[800 0 1080 1080];
    plot(500,-500,'k');
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    axis equal
    xlim([0,9000]);
    ylim([-9000,0]);
    if(playAlone == 0)
        MplayerForPlot = M.player;
    else
        MplayerForPlot = sside;
    end
    for i=1:8
        for j=1:8
            xri = 500+1000*(j-1);
            yri = -500-i*1000;
            xri = abs(8000*MplayerForPlot - xri);
            yri = -abs(10000*MplayerForPlot - abs(yri));
            xr = [xri, xri, xri+1000, xri+1000];
            yr = [yri, yri+1000, yri+1000, yri];
            if(mod(i+j,2)==0)
                fill(xr, yr, whiteBackgroundColor, 'FaceAlpha', ...
                    1.0, 'EdgeColor', 'none');
            else
                fill(xr, yr, blackBackgroundColor, 'FaceAlpha', ...
                    1.0, 'EdgeColor', 'none');
            end
            if(A(i,j)>0.5)
                if(A(i,j)==1||A(i,j)==11)
                    plotPawn(xri, yri, A(i,j)>10);
                elseif(A(i,j)==2||A(i,j)==12)
                    plotBishop(xri, yri, A(i,j)>10);
                elseif(A(i,j)==3||A(i,j)==13)
                    plotkNight(xri, yri, A(i,j)>10);
                elseif(A(i,j)==4||A(i,j)==14)
                    plotRook(xri, yri, A(i,j)>10);
                elseif(A(i,j)==5||A(i,j)==15)
                    plotQueen(xri, yri, A(i,j)>10);
                elseif(A(i,j)==6||A(i,j)==16)
                    plotKing(xri, yri, A(i,j)>10);
                end
            end
            
        end
    end
    plot([500,8500,8500,500,500], [-500, -500, -8500, -8500, -500], ...
        'k-', 'linewidth', 3);
    for i=1:7
        plot([500+1000*i,500+1000*i],[-500,-8500], 'k', 'linewidth', 1.2);
        plot([500,8500],[-1000*i-500,-1000*i-500], 'k', 'linewidth', 1.2);
    end
    aa=["a","b","c","d","e","f","g","h"];
    for i=1:8
        text(1000*i-60, -8700, aa(abs(MplayerForPlot*9-i)));
    end
    for i=1:8
        text(250,-1000*i, string(abs(MplayerForPlot*9-(9-i))));
    end
    titlestring = "";
    if(checkStatus.out_check == 1)
        i = checkStatus.out_zi;
        j = checkStatus.out_zj;
        xri = 500+1000*(j-1);
        yri = -500-i*1000;
        xri = abs(8000*MplayerForPlot - xri);
        yri = -abs(10000*MplayerForPlot - abs(yri));
        xr = [xri, xri, xri+1000, xri+1000, xri];
        yr = [yri, yri+1000, yri+1000, yri, yri];
        plot(xr,yr,'r','linewidth',1.3);
        titlestring = "Check! ";

    end
    if(checkmateStatus == 1)
        if(M.turn == 0)
            titlestring = "Checkmate! Black wins. ";
        else
            titlestring = "Checkmate! White wins. ";
        end
    end    
    if(M.staleRepetitions==1)
        titlestring = "Stalemate (3 repetitions). ";
    end
    if(M.staleMoves==1)
        titlestring = ...
            "Stalemate (75 moves with no pawn movements or captures). ";
    end
    if(staleStatus == 1)
        titlestring = "Stalemate. ";
    end
    if(checkmateStatus == 0 && staleStatus == 0 && ...
            M.staleRepetitions==0 && M.staleMoves==0)
        if(M.turn==0)
            titlestring = strcat(titlestring, "White to move.");
        else
            titlestring = strcat(titlestring, "Black to move.");
        end
        if(M.player==M.turn && playAlone==0)
            titlestring = strcat(titlestring, " It's your turn!");
        else
            titlestring = strcat(titlestring, " Stortfish is thinking...");
        end
    end
    if(mvv>0.5)
    titlestring = strcat("Move: ", string(mvv),". ",titlestring);
    end
    memory = M.lastMove;
    if(memory(1)>0.5)
        i = memory(1);
        j = memory(2);
        xri = 500+1000*(j-1);
        yri = -500-i*1000;
        xri = abs(8000*MplayerForPlot - xri);
        yri = -abs(10000*MplayerForPlot - abs(yri));
        xr = [xri, xri, xri+1000, xri+1000, xri];
        yr = [yri, yri+1000, yri+1000, yri, yri];
        plot(xr,yr,'g--','linewidth',1);
        i = memory(3);
        j = memory(4);
        xri = 500+1000*(j-1);
        yri = -500-i*1000;
        xri = abs(8000*MplayerForPlot - xri);
        yri = -abs(10000*MplayerForPlot - abs(yri));
        xr = [xri, xri, xri+1000, xri+1000, xri];
        yr = [yri, yri+1000, yri+1000, yri, yri];
        plot(xr,yr,'g-','linewidth',1.3);
        
    end
    %title(titlestring);
    if(prev==1)
    text(520,-250, titlestring, 'FontSize',12)
    end

end

% Plot a pawn
function [] = plotPawn(x0, y0, colore)
    global whitePiecesColor
    global blackPiecesColor
    xr = x0+1000*[.24,.28,.4,.35,.36,.43,.4,.46,.54,.6,.57,...
        .64,.65,.6,.72,.76,.24];
    yr = y0+1000*[.12,.31,.43,.51,.59,.65,.73,.8,.8,.73,.65,...
        .59,.51,.43,.31,.12,.12];
    plot(xr,yr,'k', 'linewidth', 1.2);
    if(colore==0)
        fill(xr, yr, whitePiecesColor, 'FaceAlpha', ...
            1.0, 'EdgeColor', 'none');
    else
        fill(xr, yr, blackPiecesColor, 'FaceAlpha', ...
            1.0, 'EdgeColor', 'none');
    end
end

% Plot a bishop
function [] = plotBishop(x0, y0, colore)
    global whitePiecesColor
    global blackPiecesColor
    xr = x0+1000*[0.13,0.2,0.48,0.33,0.33,0.37,0.3,0.3,0.5,...
        0.7,0.7,0.63,0.67,0.67,0.52,0.8,0.87,0.13];
    yr = y0+1000*[0.13,0.2,0.22,0.27,0.36,0.41,0.48,0.6,0.83,...
        0.6,0.48,0.41,0.36,0.27,0.22,0.2,0.13,0.13];
    plot(xr,yr,'k', 'linewidth', 1.2);
    if(colore==0)
        fill(xr, yr, whitePiecesColor, 'FaceAlpha', ...
            1.0, 'EdgeColor', 'none');
    else
        fill(xr, yr, blackPiecesColor, 'FaceAlpha', ...
            1.0, 'EdgeColor', 'none');    
    end
end

% Plot a knight
function [] = plotkNight(x0, y0, colore)
    global whitePiecesColor
    global blackPiecesColor
    xr = x0+1000*[0.32,0.36,0.49,0.53,0.27,0.17,0.13,0.3,0.3,...
        0.38,0.48,0.48,0.65,0.78,0.83,0.84,0.32];
    yr = y0+1000*[0.12,0.24,0.38,0.56,0.32,0.32,0.4,0.73,0.83,...
        0.77,0.83,0.78,0.72,0.59,0.37,0.12,0.12];
    plot(xr,yr,'k', 'linewidth', 1.2);
    if(colore==0)
        fill(xr, yr, whitePiecesColor, 'FaceAlpha', ...
            1.0, 'EdgeColor', 'none');
    else
        fill(xr, yr, blackPiecesColor, 'FaceAlpha', ...
            1.0, 'EdgeColor', 'none');
    end
end

% Plot a rook
function [] = plotRook(x0, y0, colore)
    global whitePiecesColor
    global blackPiecesColor
    xr = x0+1000*[0.2,0.2,0.27,0.27,0.31,0.31,.24,.24,.33,.33,...
        .44,.44,.55,.55,.67,.67,.76,.76,.69,.69,.73,.73,.8,.8,.2];
    yr = y0+1000*[0.12,.2,.2,.27,.35,.61,.68,.79,.79,.75,.75,...
        .79,.79,.75,.75,.79,.79,.68,.61,.34,.28,.2,.2,.12,.12];
    plot(xr,yr,'k', 'linewidth', 1.2);
    if(colore==0)
        fill(xr, yr, whitePiecesColor, 'FaceAlpha', ...
            1.0, 'EdgeColor', 'none');
    else
        fill(xr, yr, blackPiecesColor, 'FaceAlpha', ...
            1.0, 'EdgeColor', 'none');
    end
end

% Plot a queen
function [] = plotQueen(x0, y0, colore)
    global whitePiecesColor
    global blackPiecesColor
    xr = x0+1000*[.21,.27,.19,.15,.31,.31,.43,.50,.57,.69,.69,...
        .85,.81,.73,.79,.5,.21];
    yr = y0+1000*[.15,.38,.51,.76,.53,.84,.54,.86,.54,.84,.53,...
        .76,.51,.38,.15,.12,.15];
    plot(xr,yr,'k', 'linewidth', 1.2);
    if(colore==0)
        fill(xr, yr, whitePiecesColor, 'FaceAlpha', 1.0,...
            'EdgeColor', 'none');
    else
        fill(xr, yr, blackPiecesColor, 'FaceAlpha', 1.0,...
            'EdgeColor', 'none');
    end
end

% Plot a king
function [] = plotKing(x0, y0, colore)
    global whitePiecesColor
    global blackPiecesColor
    xr = x0+1000*[.5,.25,.25,.14,.14,.26,.43,.43,.5,.5,.43,.57,...
        .5,.5,.5,.57,.57,.74,.86,.86,.75,.75,.5];
    yr = y0+1000*[.11,.17,.33,.42,.57,.64,.57,.67,.73,.81,.81,...
        .81,.81,.86,.73,.67,.57,.64,.57,.42,.33,.17,.11];
    plot(xr,yr,'k', 'linewidth', 1.2);
    if(colore==0)
        fill(xr, yr, whitePiecesColor, 'FaceAlpha', 1.0,...
            'EdgeColor', 'none');
    else
        fill(xr, yr, blackPiecesColor, 'FaceAlpha', 1.0,...
            'EdgeColor', 'none');
    end
end

% Ask the user to make a move by selecting two proper cells
function [i, j] = collectUserMove(M)
    waitforbuttonpress
    % cp = get(gcf(),'CurrentPoint');
    % dx = 71.7;
    % dy = 71.7;
    % x0 = 159.6;
    % y0 = 126.0;
    % xcp = cp(1);
    % ycp = cp(2);
    cp = get(gca(),'CurrentPoint');
    x0 = 500;
    y0 = -8500;
    x8 = 8500;
    y8 = -500;
    dx = abs(x8-x0)/8;
    dy = abs(y8-y0)/8;
    xcp = cp(1,1);
    ycp = cp(1,2);
    j = 1+floor((xcp-x0)/dx);
    i = 9-(1+floor((ycp-y0)/dy));
    xr = [500+1000*(j-1), 500+1000*(j-1), 500+1000*(j), ...
        500+1000*(j), 500+1000*(j-1)];
    yr = [-500-i*1000, -500-i*1000+1000, -500-i*1000+1000, ...
        -500-i*1000, -500-i*1000];
    plot(xr,yr,'b','linewidth',1.3);
    pause(0.1);
    if(M.player == 1)
        i = 9-i;
        j = 9-j;
    end

end

% Ask the user to choose how to promote a pawn
function [uy_out, validChoice] = askChoice(turno)
    uy_out = 5;
    validChoice = 0;
    f2 = figure(2);
    hold on
    f2.Position = 0.8*[1030 400 600 225];
    plot(500,-500,'k');
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    xlim([0,5000]);
    ylim([-2000,0]);
    plotBishop(500, -1500, turno);
    plotkNight(1500, -1500, turno);
    plotRook(2500, -1500, turno);
    plotQueen(3500, -1500, turno);
    pause(0.5);
    waitforbuttonpress;
    % cp = get(gcf(),'CurrentPoint');
    % %display(cp);
    % dx = 74.2;
    % dy = 72.8;
    % x0 = 100.4;
    % y0 = 58.8;
    % xcp = cp(1);
    % ycp = cp(2);
    cp = get(gca(),'CurrentPoint');
    x0 = 500;
    y0 = -1500;
    x8 = 4500;
    y8 = -500;
    dx = abs(x8-x0)/4;
    dy = abs(y8-y0);
    xcp = cp(1,1);
    ycp = cp(1,2);

    j = 1+floor((xcp-x0)/dx);
    i = 2-(1+floor((ycp-y0)/dy));
    if(i==1 && j>=1 && j<=4)
        validChoice = 1;
        uy_out = 1+j;
    end
    xr = [500+1000*(j-1), 500+1000*(j-1), 500+1000*(j),...
        500+1000*(j), 500+1000*(j-1)];
    yr = [-500-i*1000, -500-i*1000+1000, -500-i*1000+1000,...
        -500-i*1000, -500-i*1000];
    plot(xr,yr,'b','linewidth',1.3);
    pause(0.5);
    close(f2);

end
