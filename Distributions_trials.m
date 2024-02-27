%function [GxT_N,VxT_N,NxT_serial] = Distributions_trials(MTSBELDG,MTBl)
%
%

function [GxT_N,VxT_N,NxT_serial] = Distributions_trials(MTSBELDG,MTBl)

ddays = 6:15;
n_days = length(ddays);
n_trials = 10;                                                  %number of trials per day
n_mice = length(unique(MTBl(:,1)));

%mattrices
MxDxTxG_PL = nan(n_mice,n_days,n_trials,25);
MxDxTxG_N = nan(n_mice,n_days,n_trials,25);
MxDxT_Nseg = nan(n_mice,n_days,n_trials);
MxDxT_door_visits = zeros(n_mice,n_days,n_trials,24);
MxDxT_n_serial = zeros(n_mice,n_days,n_trials,20);
door_gaps = -12:1:12;
for mice_ii = 1:n_mice
    for ii = 1:n_days

        %number of segment by gap size, and path length by gap size
        for jj = 1:25
            for kk = 1:n_trials
                ndx = find(MTSBELDG(:,8)==door_gaps(jj) & MTSBELDG(:,2)==kk+ddays(ii)*n_trials & MTSBELDG(:,1)==mice_ii);
                MxDxTxG_N(mice_ii,ii,kk,jj) = length(ndx); %per trial
                if ~isempty(ndx)
                    MxDxTxG_PL(mice_ii,ii,kk,jj) = mean(MTSBELDG(ndx,6));
                end                
            end
        end

        for jj = 1:n_trials
            %number of segment per trial
            ndx = find(MTSBELDG(:,1)==mice_ii & MTSBELDG(:,2)==(ddays(ii)*n_trials+jj));
            MxDxT_Nseg(mice_ii,ii,jj) = length(ndx);            

            %door visits across days
            tmp = histcounts(MTSBELDG(ndx,5),-0.5:1:23.5);
            MxDxT_door_visits(mice_ii,ii,jj,:) = [tmp(14:end) tmp(1:13)];

            %n_serial across days
            ndx = find(MTBl(:,1)==mice_ii & MTBl(:,2)==(ddays(ii)*10+jj));
            MxDxT_n_serial(mice_ii,ii,jj,:) = sum(MTBl(ndx,3:end),1);
        end        
    end
end
NsegDistxT = nan(50,n_trials);
for ii = 1:n_trials
    NsegDistxT(:,ii) = histcounts(MxDxT_Nseg(:,:,ii),0.5:1:50.5);
end
NsegDistxT = NsegDistxT./repmat(sum(NsegDistxT,1),50,1)*100;

GxT_PL = squeeze(mean(MxDxTxG_PL,[1 2],'omitnan'))';
GxT_N = squeeze(mean(MxDxTxG_N,[1 2],'omitnan'))';
GxT_N = GxT_N./repmat(sum(GxT_N,1),size(GxT_N,1),1)*100;
VxT_N = squeeze(mean(MxDxT_door_visits,[1 2],'omitnan'))';
VxT_N = VxT_N./repmat(sum(VxT_N,1),size(VxT_N,1),1)*100;
MxT = squeeze(mean(MxDxT_Nseg,2,'omitnan'));
NxT_serial = squeeze(mean(MxDxT_n_serial(:,:,:,1:20),[1 2]))';
NxT_serial = NxT_serial./repmat(sum(NxT_serial,1),size(NxT_serial,1),1)*100;


if nargout == 0
    T1 = 1;
    T2 = 8:10;
    figure;
    subplot(5,3,1);imagesc(GxT_PL);ylabel('seg length');xlabel('trials');
    subplot(5,3,2);plot(-12:1:12,mean(GxT_PL(:,T1),2,'omitnan'),'b',-12:1:12,mean(GxT_PL(:,T2),2,'omitnan'),'r');ylabel('path length (cm)');xlabel('seg length');
    subplot(5,3,4);imagesc(GxT_N);ylabel('seg length');xlabel('trials');
    subplot(5,3,5);plot(-12:1:12,mean(GxT_N(:,T1),2),'b',-12:1:12,mean(GxT_N(:,T2),2),'r');ylabel('# of segments per trial');xlabel('seg length');
    subplot(5,3,7);imagesc(VxT_N);ylabel('door#');xlabel('trials');
    subplot(5,3,8);plot(-11:1:12,mean(VxT_N(:,T1),2)/sum(mean(VxT_N(:,T1),2))*100,'b',-11:1:12,mean(VxT_N(:,T2),2)/sum(mean(VxT_N(:,T2),2))*100,'r');ylabel('% of visits');xlabel('door#');
    subplot(5,3,10);imagesc(NxT_serial);xlabel('trials');ylabel('serial bout length');
    subplot(5,3,11);plot(1:20,mean(NxT_serial(:,T1),2),'b',1:20,mean(NxT_serial(:,T2),2),'r');xlabel('serial bout length');ylabel('% of serial bouts');
    subplot(5,3,13);imagesc(NsegDistxT);xlabel('trials');ylabel('trial length');
    subplot(5,3,14);plot(1:50,mean(NsegDistxT(:,T1),2),'b',1:50,mean(NsegDistxT(:,T2),2),'r');xlabel('trial length (seg#)');ylabel('# of trials');
    subplot(5,3,15);plot(1:10,MxT,'b.',1:10,mean(MxT,1),'ro-');ylabel('# of segments');xlabel('trials');
end



