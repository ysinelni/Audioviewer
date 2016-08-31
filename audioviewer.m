function varargout = audioviewer(varargin)
% AUDIOVIEWERYDS M-file for audiovieweryds.fig
%      AUDIOVIEWERYDS, by itself, creates a new AUDIOVIEWERYDS or raises the existing
%      singleton*.
%
%      H = AUDIOVIEWERYDS returns the handle to a new AUDIOVIEWERYDS or the handle to
%      the existing singleton*.
%
%      AUDIOVIEWERYDS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUDIOVIEWERYDS.M with the given input arguments.
%
%      AUDIOVIEWERYDS('Property','Value',...) creates a new AUDIOVIEWERYDS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before audioviewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to audioviewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help audiovieweryds

% Last Modified by GUIDE v2.5 10-Aug-2016 11:56:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @audioviewer_OpeningFcn, ...
    'gui_OutputFcn',  @audioviewer_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before audiovieweryds is made visible.
function audioviewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to audiovieweryds (see VARARGIN)

% Choose default command line output for audiovieweryds
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes audiovieweryds wait for user response (see UIRESUME)
% uiwait(handles.AudioViewerYDS);
setappdata(0,'AudioViewerYDS',hObject);
Cfg = get_audio_viewer_default();
setappdata(hObject,'Cfg',Cfg);
setappdata(hObject,'Data',[]);
setappdata(hObject,'CurData',[]);
set(handles.InputPath,'string',Cfg.folder)
addpath(pwd);
try, addpath('../mmread'); catch, disp('No mmread detected'); end

%Create tab group
tmp=version; matver = str2num(tmp(1));
if matver<8,
    handles.tgroup = uitabgroup('v0','Parent', hObject,'TabLocation', 'top');
else
    handles.tgroup = uitabgroup('Parent', hObject,'TabLocation', 'top');
end
pos = get(handles.P1,'position');
uni = get(handles.P1,'units');
set(handles.tgroup,'units',uni)
set(handles.tgroup,'position',pos);
% set(handles.tgroup,'BackgroundColor',[0 0 0]);
tgrpos = get(handles.tgroup,'position');
if matver<8,
    handles.tab1 = uitab('v0','Parent', handles.tgroup, 'Title', ' Data ');
    handles.tab2 = uitab('v0','Parent', handles.tgroup, 'Title', ' Control ');
    handles.tab3 = uitab('v0','Parent', handles.tgroup, 'Title', ' Spectrum ');
    handles.tab4 = uitab('v0','Parent', handles.tgroup, 'Title', ' Run ');
    handles.tab5 = uitab('v0','Parent', handles.tgroup, 'Title', ' Org ');
    handles.tab6 = uitab('v0','Parent', handles.tgroup, 'Title', ' Test ');
else
    handles.tab1 = uitab('Parent', handles.tgroup, 'Title', ' Data ');
    handles.tab2 = uitab('Parent', handles.tgroup, 'Title', ' Control ');
    handles.tab3 = uitab('Parent', handles.tgroup, 'Title', ' Spectrum ');
    handles.tab4 = uitab('Parent', handles.tgroup, 'Title', ' Run ');
    handles.tab5 = uitab('Parent', handles.tgroup, 'Title', ' Org ');
    handles.tab6 = uitab('Parent', handles.tgroup, 'Title', ' Test ');
end
%Place panels into each tab
set(handles.P1,'Parent',handles.tab1)
set(handles.P2,'Parent',handles.tab2)
set(handles.P3,'Parent',handles.tab3)
set(handles.P4,'Parent',handles.tab4)
set(handles.P5,'Parent',handles.tab5)
set(handles.P6,'Parent',handles.tab6)
%Reposition each panel to same location as panel 1
corpos = [tgrpos(1) tgrpos(2) 0 0];
set(handles.P1,'position',get(handles.P1,'position')-corpos);
set(handles.P2,'position',get(handles.P1,'position'));
set(handles.P3,'position',get(handles.P1,'position'));
set(handles.P4,'position',get(handles.P1,'position'));
set(handles.P5,'position',get(handles.P1,'position'));
set(handles.P6,'position',get(handles.P1,'position'));

% set default parameters
set(0, 'DefaultFigureColor', 'white',...
    'DefaultAxesColor', 'white',...
    'DefaultLineLineStyle', '-',...
    'DefaultLineLineWidth', 1);

% set normalized
set( hObject, 'Units', 'Normalized' )
set( hObject, 'Resize', 'on' )
h = findobj( hObject, '-property', 'Units' );
set( h, 'Units', 'Normalized' )

function Cfg = get_audio_viewer_default()
Cfg.folder = '..\Data\';
Cfg.pressure_calibration_factor = 1;
Cfg.filter_frequencies = [];



% --- Outputs from this function are returned to the command line.
function varargout = audioviewer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1_select_data.
function pushbutton1_select_data_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1_select_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'AudioViewerYDS');
Cfg = getappdata(CurFig,'Cfg');

FilterSpec = { '*.wav'; '*.mp3' };

[InputFile, InputPath, tmp] = uigetfile(FilterSpec,'Select audio files',Cfg.folder,'MultiSelect','on');
if tmp==0, % cancel
    set(hObject,'Value',0);
    return
end
Cfg.folder = InputPath;

set(handles.InputPath,'string',Cfg.folder),drawnow;
setappdata(CurFig,'Cfg',Cfg);

if iscell(InputFile)
    BatchNo = length(InputFile);
else
    BatchNo = 1; tmp = InputFile; clear InputFile; InputFile{1} = tmp;
end

set(handles.listbox1_selected_data,'String',InputFile);
set(handles.listbox1_selected_data,'Value',1);
drawnow
listbox1_selected_data_Callback(handles.listbox1_selected_data, eventdata, handles)

% --- Executes on selection change in listbox1_selected_data.
function listbox1_selected_data_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1_selected_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'AudioViewerYDS');
Cfg = getappdata(CurFig,'Cfg');

load_data_from_file(handles)
update_screen(handles)

%%
% function show this data
function show_this_data(Data, Cfg, show_all_channels)

FigureTag = 'AudioViewer'; ViewFig=findobj('Tag',FigureTag);
if isempty(ViewFig), ViewFig=figure('Tag',FigureTag,'Name',FigureTag,'NumberTitle','off','Visible','on'); 
else, figure(ViewFig(1)); end;

Fs = Data.Fs;
NFFT = Cfg.NFFT;
chan = Cfg.chan;
PCF = Cfg.pressure_calibration_factor/20e-6;

switch Cfg.wts
    case 1 % signal
        
        if PCF>4194 & PCF<41943040, cylab = 'Pressure (Pa)'; ScaleData = PCF*20e-6;
        else, cylab = 'Pressure (a.u.)'; ScaleData = PCF*20e-6; end
        if PCF==1, cylab = 'Pressure (native units)'; ScaleData = 1; end
        
        if show_all_channels
            plot(Data.t,Data.y*ScaleData), xlabel('Time (s)'), 
        else
            plot(Data.t,Data.y(:,chan)*ScaleData), xlabel('Time (s)'), 
        end
        ylabel(cylab)
        title(sprintf('Acoustic Pressure (calibration factor = %5.2e)',PCF*20e-6))
        %         if ~isempty(Data.d)
        %             ax1 = gca;
        %             sxt = get(ax1,'xticklabel');
        %             xt = str2num(sxt);
        %             xd = interp1(Data.t,Data.d,xt);
        %             for k=1:length(xt)
        %                 if ~isnan(xd(k)), lxd(k,:) = sprintf('%s/%5.1f',sxt(k,:),xd(k)); else, lxd(k,:) = sprintf('%s/     ',sxt(k,:));  end
        %             end
        %             set(ax1,'xticklabel',lxd)
        %             xlabel('Time (s) / Distance (meters)')
        % %             ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','xticklabel',lxd)
        %
        %         end
    case 2 % Fourier
        
        f = Fs/2*linspace(0,1,NFFT/2+1);
        L = length(Data.y(:,chan));
        Y = fft(Data.y(:,chan)*PCF,NFFT)/L;
        %         NumP = floor(Cfg.PKWIDTH_S*Cfg.NFFT/SR/2);
        NumP = 1+floor(400*Cfg.NFFT/Data.Fs/2);
        B = SensitiveNonlinearIterativePeakClippingAlgorithm(abs(Y(1:NFFT/2+1)),[NumP:2:4*NumP]);
        
        plot(f/1e3,20*log10(abs(Y(1:NFFT/2+1))),f/1e3,20*log10(B)), grid
%         plot(f/1e3,abs(Y(1:NFFT/2+1)),f/1e3,B)
        xlim([Cfg.f1, Cfg.f2]), 
        title('Single-Sided Fast Fourier Transform')
        xlabel('Frequency (kHz)')
        ylabel('Amplitude (dB)')
        legend('Signal','Background')
       
    case 3 % Welch
        Hs=spectrum.welch; 
        Hs.SegmentLength=NFFT;
        psd(Hs,Data.y(:,chan)*PCF,'NFFT',NFFT,'Fs',Fs),
        xlim([Cfg.f1, Cfg.f2])
%         ylabel('dB re. 20\muPa/\surd{Hz}')
        ylabel('dB')
        
    case 4 % Periodogram
        Hs=spectrum.periodogram;
        psd(Hs,Data.y(:,chan)*PCF,'Fs',Fs),
        xlim([Cfg.f1, Cfg.f2])
        %         periodogram(Data.y(:,1),[],'onesided',512,Fs);
        
    case {5, 6, 7, 8, 9} % spectrograms

        switch Cfg.spectrogram_type
            case 1, spectrotype = 'psd'; st_label = 'Power Spectral Density (dB/Hz)';
                Hs=spectrum.welch;
                Hs.SegmentLength=NFFT;
                psdwelch = psd(Hs,Data.y(:,chan)*PCF,'NFFT',NFFT,'Fs',Fs);
            case 2, spectrotype = 'power'; st_label = 'Power Spectrum (dB)';
        end
        Nv = fix(2/3*NFFT);
        [S,F,T,P] = spectrogram(Data.y(:,chan)*PCF, hamming(NFFT), Nv, NFFT, Fs,spectrotype);
        
        [NF NT]=size(P);
        NumP = 1+floor(200*Cfg.NFFT/Data.Fs/2);
        for k=1:NT,
            B(:,k) = SensitiveNonlinearIterativePeakClippingAlgorithm(P(:,k),[NumP:2:4*NumP]);
        end

        %         B = repmat(mean(B,2),1,NT);
        
        if Cfg.wts<8,
            switch Cfg.wts
                case 5, X = P;
                case 6, X = B;
                case 7, X = P-B+eps;
            end
            
            surf(Data.t(1)+T,F/1e3,10*log10(abs(X)),'EdgeColor','none');
            axis xy; axis tight; colormap(jet); view(0,90);
            ylim([Cfg.f1, Cfg.f2])
            if ~Cfg.dbauto, set(gca,'clim',[Cfg.db1 Cfg.db2]); end
            colorbar
            title(st_label)
            xlabel('Time (s)'), ylabel('Frequency (kHz)'),
            
%             figure(22), spectrogram(Data.y(:,chan)*PCF, hamming(NFFT), Nv, NFFT, Fs,'power');
%             view(90,270), if ~Cfg.dbauto, set(gca,'clim',[Cfg.db1 Cfg.db2]); end
%             
%             figure(34), plot(F/1e3,min(10*log10(abs(X)),[],2),...
%                 F/1e3,mean(10*log10(abs(X)),2),...
%                 F/1e3,max(10*log10(abs(X)),[],2),...
%                 psdwelch.Frequencies/1e3, 10*log10(psdwelch.Data)), grid
%             xlabel('Frequency (kHz)'), ylabel(st_label), legend('min','average','max','welch')

        else
            switch Cfg.wts
                case 8, 
                    % phase diff along frequency axis
                    NDIFF =1;
                    surf(Data.t(1)+T,F(1:end-NDIFF)/1e3,diff(unwrap(angle(S),[],1),NDIFF,1),'EdgeColor','none');
                    title('aka Group delay')
                case 9, 
                    % phase diff along time axis
                    NDIFF = 1;
                    surf(Data.t(1)+T(1:end-NDIFF),F/1e3,diff(unwrap(angle(S),[],2),NDIFF,2),'EdgeColor','none');
                    title('aka Instantaneous frequency')
            end
            axis xy; axis tight; colormap(jet); view(0,90);
            ylim([Cfg.f1, Cfg.f2])
            colorbar
            xlabel('Time (s)'), ylabel('Frequency (kHz)'),
            
            switch Cfg.wts
                case 8,
                    display_hough_transform(diff(unwrap(angle(S),[],1),NDIFF,1))
                case 9,
                    display_hough_transform(diff(unwrap(angle(S),[],2),NDIFF,2))
            end
        end
        
    case 10 % Distance
        try
%         if ~isempty(Data.d)
            plot(Data.t,Data.d,'linewidth',2), xlabel('Time (s)'), ylabel('Distance (m)'), grid
%         end
        catch
            cla; title('Distance')
            disp('Error: no distance provided')
        end
        
end

set(gca,'fontsize',32);
title('')


function display_hough_transform(I)
% Convert to intensity.
% I  = rgb2gray(RGB);

figure(33)
se = strel('line',20,0);
ac = imclose(I,se);
figure(33), imagesc(ac), colorbar
set(gca,'ydir','normal')
colormap(bone);


% % Extract edges.
% BW = edge(I,'canny',0.7);
% [H,T,R] = hough(BW,'RhoResolution',0.5,'Theta',-90:0.5:89.5);
% 
% figure(33)
% % Display the original image.
% subplot(2,1,1);
% imshow(BW);
% title('aka Instantaneous frequency')
% 
% % Display the Hough matrix.
% subplot(2,1,2);
% imshow(imadjust(mat2gray(H)),'XData',T,'YData',R,...
%       'InitialMagnification','fit');
% title('Hough Transform');
% xlabel('\theta'), ylabel('\rho');
% axis on, axis normal, hold on;
% colormap(hot);


% --- Executes during object creation, after setting all properties.
function listbox1_selected_data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1_selected_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function data_duration_Callback(hObject, eventdata, handles)
% hObject    handle to data_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of data_duration as text
%        str2double(get(hObject,'String')) returns contents of data_duration as a double


% --- Executes during object creation, after setting all properties.
function data_duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function data_sampling_frequency_Callback(hObject, eventdata, handles)
% hObject    handle to data_sampling_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of data_sampling_frequency as text
%        str2double(get(hObject,'String')) returns contents of data_sampling_frequency as a double


% --- Executes during object creation, after setting all properties.
function data_sampling_frequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_sampling_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function data_channels_Callback(hObject, eventdata, handles)
% hObject    handle to data_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of data_channels as text
%        str2double(get(hObject,'String')) returns contents of data_channels as a double


% --- Executes during object creation, after setting all properties.
function data_channels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2_what_to_show.
function listbox2_what_to_show_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2_what_to_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2_what_to_show contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2_what_to_show

w = get(handles.listbox2_what_to_show,'value');
if w>=5 & w<=9
    set(handles.spectrogram_type,'Enable','on')
else
    set(handles.spectrogram_type,'Enable','off')
end    
try
    update_screen(handles)
catch
    disp('Error updating screen - load data first')
end

%%
% Update screen
function update_screen(handles)
CurFig = getappdata(0,'AudioViewerYDS');
Cfg = getappdata(CurFig,'Cfg');
CurData = getappdata(CurFig,'CurData');

try, Cfg.wtsold = Cfg.wts; catch, Cfg.wtsold=0; end
Cfg.wts = get(handles.listbox2_what_to_show,'Value');
Cfg.chan = get(handles.show_channel,'value');
Cfg.NFFT = round(str2double(get(handles.NFFT,'String')));

Cfg.f1 = str2double(get(handles.show_start_freq,'string'));
Cfg.f2 = str2double(get(handles.show_stop_freq,'string'));

Cfg.db1 = str2double(get(handles.spectr_db_min,'string'));
Cfg.db2 = str2double(get(handles.spectr_db_max,'string'));
Cfg.dbauto = get(handles.spectr_autoscale,'value');

Cfg.pressure_calibration_factor = str2double(get(handles.pressure_calibration_factor,'string'));
Cfg.spectrogram_type = get(handles.spectrogram_type,'value');

spl_val = calculate_spl_values(CurData.y(:,Cfg.chan),CurData.Fs,Cfg.pressure_calibration_factor);

show_all_channels = get(handles.show_all_channels,'value');

set(handles.spl_dba,'String',sprintf('%5.1f',spl_val(1)));

setappdata(CurFig,'Cfg',Cfg);
% CurData.d = [];
% cleardistfigflag = 1;
% FigureTag = 'Debug distances'; hfig=findobj('Tag',FigureTag);
% if isempty(hfig), hfig=figure('Tag',FigureTag,'Name',FigureTag,'NumberTitle','off'); end;
% 
% if cleardistfigflag, close(hfig); end

sfigure(CurFig)
show_this_data(CurData, Cfg, show_all_channels)

%%
% Load data from file
function load_data_from_file(handles)
CurFig = getappdata(0,'AudioViewerYDS');
Cfg = getappdata(CurFig,'Cfg');
Data = getappdata(CurFig,'Data');

contents = cellstr(get(handles.listbox1_selected_data,'String'));
selected_value = get(handles.listbox1_selected_data,'Value');
selected = contents{selected_value};
FileName = fullfile(Cfg.folder,selected);
Cfg.filename = selected;
setappdata(CurFig,'Cfg',Cfg);
tmp=version; matver = str2num(tmp(1));

FileType = -1;
if strcmp(lower(FileName(end-2:end)),'wav'), FileType = 1; end
if strcmp(lower(FileName(end-2:end)),'mp3'), FileType = 2; end

down_sample_factor = str2double(get(handles.down_sample_factor,'string'));
if down_sample_factor <1 | floor(down_sample_factor)~=ceil(down_sample_factor), 
    disp('Error: down sample factor must be integer biger than 1'), return; end

switch FileType
    case 1
        if matver<8,
            [m d] = wavfinfo(FileName);
        else
            info = audioinfo(FileName); m=info.NumChannels;
        end
        if isempty(m),
            disp('Error reding data file');
            try_to_fix_wav_file(FileName);
            return,
        end
        if matver<8,
            sizeinfo = wavread(FileName,'size');
            [y, fs, nbits, opts] = wavread(FileName,1);
        else
            sizeinfo = [info.TotalSamples, info.NumChannels];
            [y, fs] = audioread(FileName,[1 1]);
            nbits = info.BitsPerSample;
        end
    case 2
        [v,a]=mmread(FileName,0,[],0);
        fs = a.rate;
        nbits = a.bits;
        sizeinfo = [a.totalDuration*fs, a.nrChannels];
    otherwise
        disp('Error: unrecognized file type.')
        return
end

set(handles.data_sampling_frequency,'string',sprintf('%5.2f',fs/1e3))


freq_autoscale = get(handles.freq_autoscale,'value');

if freq_autoscale
    set(handles.show_stop_freq,'string',sprintf('%6.3f',fs/2e3))
    set(handles.show_start_freq,'string',sprintf('%6.3f',0))
end

data_duration = sizeinfo(1) / fs;
set(handles.data_duration,'string',sprintf('%5.2f',data_duration))

data_channels = sizeinfo(2);
set(handles.data_channels,'string',sprintf('%d',data_channels))
for k=1:data_channels, chs{k} = sprintf('%02d',k); end
set(handles.show_channel,'string',chs)
chan = get(handles.show_channel,'value');
if chan>data_channels,
    set(handles.show_channel,'value',1)
end

fpos = get(CurFig,'position');
h = msgbox(sprintf('Loading %dx%d %d bits audio data',sizeinfo(1),sizeinfo(2),nbits),'Wait','warn');
set(h,'units',get(CurFig,'units'))
mpos = get(h,'position');
set(h,'position',[fpos(1)+(fpos(3)-mpos(3))/2 fpos(2)+(fpos(4)-mpos(4))/2 mpos(3), mpos(4)])
ob = findobj(h,'tag','OKButton'); delete(ob);
mb = findobj(h,'tag','MessageBox');
drawnow

switch FileType
    case 1
        if matver<8,
            [Y, fs] = wavread(FileName, 'double');
        else
            [Y, fs] = audioread(FileName, 'double');
%             Y = Y*2^(nbits-1); % 8/22 no longer needed after calibration
%             factor has been changed to Pa/V definition
%             Y=Y/256; % to make it compatible with old wavread, only
%             applicable if reading native into 32 bits
        end
    case 2
        [v,a]=mmread(FileName,0,[],0);
        Y = a.data;
        clear a v;
end

if down_sample_factor==1
    Data.Y = Y;
    Data.Fs = fs;
else
    Data.Y = Y(1:down_sample_factor:end,:);
    Data.Fs = fs/down_sample_factor;
    set(handles.data_sampling_frequency,'string',sprintf('%5.2f',Data.Fs/1e3))
end


setappdata(CurFig,'Data',Data)
drawnow;

t1 = str2double(get(handles.show_start_time,'string'));
t2 = t1+str2double(get(handles.show_window_time,'string'));

n1 = floor(t1 * Data.Fs) + 1;
n2 = floor(t2 * Data.Fs);

[NT,NC] = size(Data.Y);
if n1 > NT, n1 = NT; end
if n2 > NT, n2 = NT; end


%%
% see if organizer was loaded
try
    chan_num = get(handles.show_channel,'value');
    Cfg.pressure_calibration_factor = Cfg.org(selected_value).calibration(chan_num);
    
    disp(sprintf('Calibration for chanel %d is set to %5.2f from organizer',chan_num, Cfg.pressure_calibration_factor));
    set(handles.pressure_calibration_factor,'string',sprintf('%5.2f',Cfg.pressure_calibration_factor)),
    setappdata(CurFig,'Cfg',Cfg), drawnow;
catch
    disp(sprintf('Calibration was not loaded from organizer: assuming %5.2f',Cfg.pressure_calibration_factor));
end

PCF = Cfg.pressure_calibration_factor;
for k=1:NC
    y1 = double(max(Data.Y(n1:n2,k)));
    y2 = double(min(Data.Y(n1:n2,k)));
    
    disp(sprintf('Raw data channel %d [%4.1f-%4.1f] s. max min: %4.2e %4.2e',k,t1,t2,y1,y2))
    disp(sprintf('                        corrected max min: %4.2f %4.2f Pa',y1*PCF,y2*PCF))
end
update_current_data(t1,t2)
drawnow;

try, close(h); catch, end


%%
% try to fix wav file
function try_to_fix_wav_file(FullFileName)
fprintf('Attempting to fix chunk size in the wave file: ')
try
    wavchunksizefix( FullFileName );
    fprintf('Done\n')
catch
    fprintf('Could not fix\n')
end


%%
% Upload data from memory
function update_current_data(t1,t2)
% disp('function update_current_data')
% tic
CurFig = getappdata(0,'AudioViewerYDS');
Data = getappdata(CurFig,'Data');
CurData = getappdata(CurFig,'CurData');
Cfg = getappdata(CurFig,'Cfg');


% toc
n1 = floor(t1 * Data.Fs) + 1;
n2 = floor(t2 * Data.Fs);

[NT,NC] = size(Data.Y);

if n1 > NT, n1 = NT; end
if n2 > NT, n2 = NT; end

CurData.Fs = Data.Fs;
CurData.t = (n1:n2)/Data.Fs;
CurData.y = double(Data.Y(n1:n2,:));

if ~isempty(Cfg.filter_frequencies)
    CurData.y = FilterSingal(CurData.y,Cfg,CurData.Fs);
end

% toc
setappdata(CurFig,'CurData',CurData), drawnow;
% toc

% --- Executes during object creation, after setting all properties.
function listbox2_what_to_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2_what_to_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3_show.
function pushbutton3_show_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

t1 = str2double(get(handles.show_start_time,'string'));
t2 = t1+str2double(get(handles.show_window_time,'string'));

update_current_data(t1,t2)
update_screen(handles)


function show_start_time_Callback(hObject, eventdata, handles)
% hObject    handle to show_start_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of show_start_time as text
%        str2double(get(hObject,'String')) returns contents of show_start_time as a double
t1 = str2double(get(handles.show_start_time,'string'));
t2 = t1+str2double(get(handles.show_window_time,'string'));

update_current_data(t1,t2)
update_screen(handles)

% --- Executes during object creation, after setting all properties.
function show_start_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to show_start_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function show_window_time_Callback(hObject, eventdata, handles)
% hObject    handle to show_window_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of show_window_time as text
%        str2double(get(hObject,'String')) returns contents of show_window_time as a double
t1 = str2double(get(handles.show_start_time,'string'));
t2 = t1+str2double(get(handles.show_window_time,'string'));

update_current_data(t1,t2)
update_screen(handles)

% --- Executes during object creation, after setting all properties.
function show_window_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to show_window_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NFFT_Callback(hObject, eventdata, handles)
% hObject    handle to NFFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NFFT as text
%        str2double(get(hObject,'String')) returns contents of NFFT as a double
update_screen(handles)

% --- Executes during object creation, after setting all properties.
function NFFT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NFFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in show_channel.
function show_channel_Callback(hObject, eventdata, handles)
% hObject    handle to show_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns show_channel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from show_channel
update_screen(handles)

% --- Executes during object creation, after setting all properties.
function show_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to show_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function show_start_freq_Callback(hObject, eventdata, handles)
% hObject    handle to show_start_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of show_start_freq as text
%        str2double(get(hObject,'String')) returns contents of show_start_freq as a double
update_screen(handles)

% --- Executes during object creation, after setting all properties.
function show_start_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to show_start_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function show_stop_freq_Callback(hObject, eventdata, handles)
% hObject    handle to show_stop_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of show_stop_freq as text
%        str2double(get(hObject,'String')) returns contents of show_stop_freq as a double
update_screen(handles)

% --- Executes during object creation, after setting all properties.
function show_stop_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to show_stop_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4_play_signal.
function pushbutton4_play_signal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4_play_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'AudioViewerYDS');
Cfg = getappdata(CurFig,'Cfg');
CurData = getappdata(CurFig,'CurData');

y = CurData.y(:,Cfg.chan);
Fs = CurData.Fs;

scff = str2double(get(handles.scale_audio_frequency_factor,'string'));

update_screen(handles),

if Cfg.wts~=5 & Cfg.wts~=1,
    disp('Player is only available for signal and spectrogram')
    return
end
% hold(handles.mainaxis,'on')
hold('on')

set(handles.pushbutton4_play_signal,'enable','off')

zlimits = get( 'ZLim'); % get the y-axis limits
ylimits = get( 'YLim'); % get the y-axis limits
% plotdata = [ylimits(1):0.1:ylimits(2)];
hline = plot(CurData.t(1) + repmat(0, size(ylimits)), ylimits,'g',...
    'tag','plotmarker','ZData',repmat(zlimits(2)+1, size(zlimits))); % plot the marker

%% instantiate the audioplayer object
player = audioplayer(y/max(y), Fs*scff);


%% set conversion and axis
plotdata.ylimits = ylimits;
plotdata.zlimits = zlimits;
% plotdata.axis = handles.mainaxis;
plotdata.t1 = CurData.t(1);
plotdata.s2t = abs(CurData.t(end)-CurData.t(1))/player.TotalSamples;


%% setup the timer for the audioplayer object
player.TimerFcn = {@plotMarker, player, CurFig, plotdata}; % timer callback function (defined below)
player.TimerPeriod = 0.1; % period of the timer in seconds

%% start playing the audio
% this will move the marker over the audio plot at intervals of 0.01 s
play(player);

while isplaying(player)
    pause(1)
end
try, hMarker = findobj(CurFig, 'tag','plotmarker');
    delete(hMarker); catch, end
% hold(handles.mainaxis,'off')
hold('off')
%disp('done playing')
set(handles.pushbutton4_play_signal,'enable','on')
% wavwrite(y/max(y),Fs,'audioviewer.wav')
audiowrite('audioviewer.wav',y,Fs)


function spectr_db_min_Callback(hObject, eventdata, handles)
% hObject    handle to spectr_db_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spectr_db_min as text
%        str2double(get(hObject,'String')) returns contents of spectr_db_min as a double
update_screen(handles)

% --- Executes during object creation, after setting all properties.
function spectr_db_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spectr_db_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function spectr_db_max_Callback(hObject, eventdata, handles)
% hObject    handle to spectr_db_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spectr_db_max as text
%        str2double(get(hObject,'String')) returns contents of spectr_db_max as a double
update_screen(handles)

% --- Executes during object creation, after setting all properties.
function spectr_db_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spectr_db_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in spectr_autoscale.
function spectr_autoscale_Callback(hObject, eventdata, handles)
% hObject    handle to spectr_autoscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of spectr_autoscale
update_screen(handles)

% --- Executes on button press in freq_autoscale.
function freq_autoscale_Callback(hObject, eventdata, handles)
% hObject    handle to freq_autoscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of freq_autoscale



function scale_audio_frequency_factor_Callback(hObject, eventdata, handles)
% hObject    handle to scale_audio_frequency_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scale_audio_frequency_factor as text
%        str2double(get(hObject,'String')) returns contents of scale_audio_frequency_factor as a double


% --- Executes during object creation, after setting all properties.
function scale_audio_frequency_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scale_audio_frequency_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% ------------------------------------------------------------------------
%% the timer callback function definition
function plotMarker(...
    obj, ...            % refers to the object that called this function (necessary parameter for all callback functions)
    eventdata, ...      % this parameter is not used but is necessary for all callback functions
    player, ...         % we pass the audioplayer object to the callback function
    figHandle, ...      % pass the figure handle also to the callback function
    plotdata)           % finally, we pass the data necessary to draw the new marker

% check if sound is playing, then only plot new marker
if strcmp(player.Running, 'on')
    %     disp('this one')
    % get the handle of current marker and delete the marker
    hMarker = findobj(figHandle,'tag','plotmarker');
    delete(hMarker);
    
    % get the currently playing sample
    x = plotdata.t1 + plotdata.s2t*player.CurrentSample;
    
    % plot the new marker
    h=plot(plotdata.axis,repmat(x, size(plotdata.ylimits)), plotdata.ylimits, 'g',...
        'tag','plotmarker','ZData',repmat(plotdata.zlimits(2)+1, size(plotdata.zlimits)));
    %     get(h)
    
    
end

%%
% Background
function vv = SensitiveNonlinearIterativePeakClippingAlgorithm(v,P)
%%
%% v - vector
%% P clipping order value or vector
%%
N = length(v);

for p=P,
    for n=1:N,
        if n-p<1,
%             vv(n) = v(n);
            vv(n) = min(v(n),v(n+p));
        elseif n+p>N
%             vv(n) = v(n);
            vv(n) = min(v(n),v(n-p));
        else
            vv(n) = min(v(n),(v(n+p)+v(n-p))/2);
        end
    end
    v = vv;
end
vv=vv';

%%
% prepare figure to be updated in the background
function sfigure(h)
if nargin>=1
    if ishandle(h), set(0, 'CurrentFigure', h);
    else h = figure(h); end %#ok<*NASGU>
else h = figure;
end

%%
% fix brocken wav file
function wavchunksizefix( filename )
d = dir(filename);
fileSize = d.bytes;
fid=fopen(filename,'r+','l');
fseek(fid,4,-1);
fwrite(fid,fileSize-8,'uint32');
fseek(fid,40,-1);
fwrite(fid,fileSize-44,'uint32');
fclose(fid);



function spl_dba_Callback(hObject, eventdata, handles)
% hObject    handle to spl_dba (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spl_dba as text
%        str2double(get(hObject,'String')) returns contents of spl_dba as a double


% --- Executes during object creation, after setting all properties.
function spl_dba_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spl_dba (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function spl_vals =  calculate_spl_values(x,Fs,pcf)
% x - signal in time domain
% Fs - sampling frequency
% pcf - pressure callibration factor
REFERENCE_SOUND_PRESSURE = 20e-6; 
NP = length(x);
NFFT = length(x);
Duration = NP/Fs;
x = x*pcf;

f = Fs*(-NFFT/2:NFFT/2-1)/NFFT;
t = (1:NP)/Fs;

X = fft(x,NFFT);
equivalent_sound_level = 20*log10(rms(x)/REFERENCE_SOUND_PRESSURE);

spl_f = 10*log10( (abs(fftshift(X))).^2/Fs/NFFT ); % double check
            
avpwr_f = trapz(f,(abs(fftshift(X))).^2)/NFFT/Fs;
avpwr_t = trapz((x.^2) )/Duration/Fs;
avpwr_t = mean(x.^2); % this should be correct per disc with A. Sedunov 7/7/16

esl_f = 20*log10(sqrt(avpwr_f)/REFERENCE_SOUND_PRESSURE);
esl_t = 20*log10(sqrt(avpwr_t)/REFERENCE_SOUND_PRESSURE);

%%
% a-weighted filter

% filter coefficients
c1 = 3.5041384e16;
c2 = 20.598997^2;
c3 = 107.65265^2;
c4 = 737.86223^2;
c5 = 12194.217^2;

% evaluate A-weighting filter
ff = f.^2;
num = c1*ff.^4;
den = ((c2+ff).^2) .* (c3+ff) .* (c4+ff) .* ((c5+ff).^2);
A = num./den;

% filtering
sa = spl_f.*A.';

% %%
% % calculate loudness - need to double check
% if f(1)<0, % this is two sided distribution
%     saval = 20*log10( (trapz( abs(diff(f(1:2)))*10.^(sa/10) ))/NFFT/Fs/2 )/2;
% else
%     saval = 20*log10( (trapz( abs(diff(f(1:2)))*10.^(sa/10) ))/NFFT/Fs/2 );
% end

[LevelAW LevelKhz] = aweighting_sedunov(x,Fs,REFERENCE_SOUND_PRESSURE);

spl_vals = [LevelAW, esl_f, esl_t];

function [LevelAW LevelKhz] = aweighting_sedunov(x,Fs,REFERENCE_SOUND_PRESSURE)
% A-weighting function realization
AweightFormula = @(f) 1.2588966*(12200^2.*f.^4)./ ...
    ((f.^2+20.6^2).*sqrt((f.^2+107.7^2).*(f.^2+737.9^2)).*(f.^2+12200^2));

% compute psd

[Pxx,W]=pwelch(x);
F = linspace(0,Fs/2,numel(W));
SLR = F(2); %spectral line resolution

%  compute a-weighting for the frequency range

AW = AweightFormula(F);

WOnly1kHz = (F<1020) & (F>980);

%  This is a-weighted PSD
%  this is what a power meter would measure
PxxAW = Pxx.*(AW(:));

% integrate the power over frequencies

% 20*log10(rms(x))
Level = db(sqrt(sum(Pxx)*W(2)))-db(REFERENCE_SOUND_PRESSURE);
LevelAW = db(sqrt(sum(PxxAW)*W(2)))-db(REFERENCE_SOUND_PRESSURE);
LevelKhz = db(sqrt(sum(PxxAW.*WOnly1kHz(:))*W(2)))-db(REFERENCE_SOUND_PRESSURE);
            
function y = rms(u)

if size(u,1) == 1,
    %     disp('rms warning: size is one'),
    flip = 1;
    u=u.'; else, flip=0; end
y = sqrt(sum(u.*conj(u))/size(u,1));
if flip, y=y'; end


function pressure_calibration_factor_Callback(hObject, eventdata, handles)
% hObject    handle to pressure_calibration_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pressure_calibration_factor as text
%        str2double(get(hObject,'String')) returns contents of pressure_calibration_factor as a double
update_screen(handles)

% --- Executes during object creation, after setting all properties.
function pressure_calibration_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pressure_calibration_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6_run.
function pushbutton6_run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'AudioViewerYDS');
Cfg = getappdata(CurFig,'Cfg');
Data = getappdata(CurFig,'Data');

Fs = Data.Fs;
NFFT = Cfg.NFFT;
chan_num = get(handles.show_channel,'value');


REFERENCE_SOUND_PRESSURE = 20e-6;
Z_AIR = 410; % Pa·s/m = kg /m^2 / s
                
tmp=version; matver = str2num(tmp(1));

FileNames = get(handles.listbox1_selected_data,'string'); 
NumberOfFiles = length(FileNames);
RunValue = get(handles.pushbutton6_run,'value');
if ~RunValue, return; end

run_type = get(handles.run_type,'value');

run_state = get(handles.run_state,'backgroundcolor');

set(handles.run_state,'backgroundcolor',[0 1 0]);

run_range_min = str2double(get(handles.run_range_min,'string'));
run_range_max = str2double(get(handles.run_range_max,'string'));
run_overlap_percent  = str2double(get(handles.run_overlap_percent,'string'));

ThisRunStartTime = tic;
% initialize statistics
Stat.Cfg=Cfg; stat_num=1;

selected_file = get(handles.listbox1_selected_data,'Value')
        
for k = selected_file:NumberOfFiles;
    set(handles.listbox1_selected_data,'Value',k), drawnow;
    load_data_from_file(handles)

    try
        CurData = getappdata(CurFig,'CurData');
        CurCalFactor = Cfg.org(k).calibration(1:length(CurData.y(1,:)));
        
        disp(sprintf('Calibration for chanel %d is set to %5.2f from organizer',chan_num, CurCalFactor(chan_num)));
        set(handles.pressure_calibration_factor,'string',sprintf('%5.2e',CurCalFactor(chan_num))), 
        Cfg.pressure_calibration_factor = CurCalFactor(chan_num);
        setappdata(CurFig,'Cfg',Cfg), drawnow;
    catch
        disp(sprintf('Calibration was not loaded from organizer: assuming %4.3e',Cfg.pressure_calibration_factor));
        CurCalFactor = ones([length(CurData.y(1,:)),1])*Cfg.pressure_calibration_factor;
    end
    
    
    Stat.CurCalFactor = CurCalFactor;

    
    data_duration = str2double(get(handles.data_duration,'string'));
    show_window_time = str2double(get(handles.show_window_time,'string'));
    increment_time = show_window_time*(100-run_overlap_percent)/100;
    
    
    tmp = char(FileNames{k}); 
    outfilemask = fullfile(Cfg.folder,[tmp(1:end-4),sprintf('_%03d',chan_num)]);
    outstatmask = fullfile(Cfg.folder,tmp(1:end-4));
    
    frame_seq = 1;
    frame_num = 1; clear Frame;
    frame_seq_start = 0;
    for show_start_time = 0:increment_time:(data_duration-show_window_time)
        
        cureltime = toc(ThisRunStartTime);
        set(handles.show_start_time,'string',sprintf('%8.4f',show_start_time))
        set(handles.run_status,'string',sprintf('%s running %s sequence %d starting at %7.4f seconds',...
            datestr(cureltime/86400,13), char(FileNames{k}), k, show_start_time)), drawnow;
        
        t1 = str2double(get(handles.show_start_time,'string'));
        t2 = t1+show_window_time;
        
        update_current_data(t1,t2)
        update_screen(handles)
        
        switch Cfg.wts
            case 1 % signal
                ylim([run_range_min run_range_max])
            case {5, 6, 7, 8, 9} % spectrogram
        end
        
        switch run_type
            case 1,% calculate statistics
                CurData = getappdata(CurFig,'CurData'); drawnow;
                Stat.window_time = show_window_time;
                Stat.start_time(stat_num) = show_start_time;
                for nncc=1:length(CurData.y(1,:))
                    Stat.avpwr_t(stat_num,nncc) = mean((CurCalFactor(nncc)*CurData.y(:,nncc)).^2);
                    Stat.esl_t(stat_num,nncc) = 20*log10(sqrt(Stat.avpwr_t(stat_num,nncc))/REFERENCE_SOUND_PRESSURE);
                    [LevelAW, LevelKhz] = aweighting_sedunov(CurCalFactor(nncc)*CurData.y(:,nncc),CurData.Fs,REFERENCE_SOUND_PRESSURE);
                    Stat.dba(stat_num,nncc) = LevelAW;
                    Stat.dba_1kHz(stat_num,nncc) = LevelKhz;
                end
                stat_num=stat_num+1;
                clear CurData;
                
            case 2, % save movies
                Frame(frame_num) = getframe(gca); drawnow;
                frame_num=frame_num+1;
                
                if frame_num>256 | (show_start_time > (data_duration-show_window_time-increment_time)) | ~RunValue,
                    
                    outfile_avi = [outfilemask,sprintf('_%03d',frame_seq),'.avi'];
                    if matver<8
                        movie2avi(Frame,outfile_avi); pause(1);
                    else
                        try
                            writerObj = VideoWriter(outfile_avi);
                            open(writerObj);
                            writeVideo(writerObj,Frame); pause(1);
                            close(writerObj);
                            clear writerObj;
                        catch
                            disp('Ups: could not write video')
                        end
                    end
                    
                    [y,Fs] = get_section_of_audio_data(frame_seq_start,show_start_time + increment_time);
                    outfile_wav = [outfilemask,sprintf('_%03d',frame_seq),'.wav'];
                    if matver<8
                        wavwrite(y/max(y),Fs,outfile_wav); pause(0.1);
                    else
                        audiowrite(outfile_wav,y,Fs); pause(0.1);
                    end
                    
                    frame_seq_start = show_start_time + increment_time;
                    
                    clear Frame y Fs;
                    frame_num = 1;
                    frame_seq = frame_seq + 1;
                    
                end
            case 3,% calculate psd and energy

                CurData = getappdata(CurFig,'CurData'); drawnow;
                Stat.Fs = Fs;
                Stat.window_time = show_window_time;
                Stat.calibration_factor = CurCalFactor(chan_num);
                Stat.start_time(stat_num) = show_start_time;
                
                Hs=spectrum.welch;
                Hs.SegmentLength=NFFT;
                PCF = CurCalFactor(chan_num)/REFERENCE_SOUND_PRESSURE;
                
                x = CurData.y(:,chan_num)*PCF;
                
                psdwelch = psd(Hs,x, 'NFFT',NFFT, 'Fs',Fs);
                
                
                
                Stat.energy(stat_num) = trapz(x.^2)/Fs;
                Stat.power(stat_num) = Stat.energy(stat_num)/show_window_time;
                
                Stat.psd(stat_num,:) = psdwelch.Data;
                Stat.f = psdwelch.Frequencies;
                
                Stat.var_sig(stat_num) = var(x);
                Stat.psd_sum_f(stat_num) = trapz(psdwelch.Data)*abs(Stat.f(2)-Stat.f(1));
                
                stat_num=stat_num+1;
                
                clear CurData;
        end
        RunValue = get(handles.pushbutton6_run,'value'); drawnow;
        if ~RunValue, break; end
    end
    
    switch run_type
        case 1, % save statistics
            outfile_stat = [outstatmask,' stat.mat'];
            save(outfile_stat, 'Stat'); pause(0.1);
            stat_num=1; clear Stat;
        case 3, % save psd
            outfile_stat = [outstatmask,' psd.mat'];
            save(outfile_stat, 'Stat'); pause(0.1);
            stat_num=1; clear Stat;
    end
    set(handles.run_status,'string',sprintf('%s completed %s sequence %d',...
            datestr(cureltime/86400,13), char(FileNames{k}))), drawnow;
    set(handles.show_start_time,'string','0')
    
    if ~RunValue, break; end
end


set(handles.run_state,'backgroundcolor',[0.3 0.75 0.93]);
set(handles.pushbutton6_run,'value',0);

%%
% return section of audio data corresponding to time window
function [y,Fs] = get_section_of_audio_data(t1,t2)
CurFig = getappdata(0,'AudioViewerYDS');
Data = getappdata(CurFig,'Data');
Cfg = getappdata(CurFig,'Cfg');

n1 = floor(t1 * Data.Fs) + 1;
n2 = floor(t2 * Data.Fs);

[NT,NC] = size(Data.Y);

if n1 > NT, n1 = NT; end
if n2 > NT, n2 = NT; end

Fs = Data.Fs;
y = double(Data.Y(n1:n2,Cfg.chan));


function run_range_min_Callback(hObject, eventdata, handles)
% hObject    handle to run_range_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of run_range_min as text
%        str2double(get(hObject,'String')) returns contents of run_range_min as a double


% --- Executes during object creation, after setting all properties.
function run_range_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to run_range_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function run_range_max_Callback(hObject, eventdata, handles)
% hObject    handle to run_range_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of run_range_max as text
%        str2double(get(hObject,'String')) returns contents of run_range_max as a double


% --- Executes during object creation, after setting all properties.
function run_range_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to run_range_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function run_overlap_percent_Callback(hObject, eventdata, handles)
% hObject    handle to run_overlap_percent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of run_overlap_percent as text
%        str2double(get(hObject,'String')) returns contents of run_overlap_percent as a double
run_overlap_percent  = str2double(get(handles.run_overlap_percent,'string'));
if run_overlap_percent<0 | run_overlap_percent >=100
    run_overlap_percent = 0;
    set(handles.run_overlap_percent,'string','0.0')
end
    
% --- Executes during object creation, after setting all properties.
function run_overlap_percent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to run_overlap_percent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in show_all_channels.
function show_all_channels_Callback(hObject, eventdata, handles)
% hObject    handle to show_all_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_all_channels
update_screen(handles)


% --- Executes on selection change in run_events.
function run_events_Callback(hObject, eventdata, handles)
% hObject    handle to run_events (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'AudioViewerYDS');
Cfg = getappdata(CurFig,'Cfg');

events = cellstr(get(handles.run_events,'String')) ;
if isempty(events), return; end
selected_event = events{get(handles.run_events,'Value')};

event_number = str2double(get(handles.event_number,'String'));

for k=1:length(Cfg.unique_events.name)
    if strcmp(char(Cfg.unique_events.name{k}),selected_event), break, end
end
set(handles.events_text,'string',sprintf('Event %d %s has occurence of %d',k, upper(char(Cfg.unique_events.name{k})),Cfg.unique_events.occurence(k)))

if event_number>Cfg.unique_events.occurence(k),
    disp(sprintf('Warning: event index %s exceeds occurence and is changed to 1',upper(char(Cfg.unique_events.name{k}))))
    set(handles.event_number,'String','1')
end

% --- Executes during object creation, after setting all properties.
function run_events_CreateFcn(hObject, eventdata, handles)
% hObject    handle to run_events (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7_read_events.
function pushbutton7_read_events_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7_read_events (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'AudioViewerYDS');
Data = getappdata(CurFig,'Data');
Cfg = getappdata(CurFig,'Cfg');

[File2read, Path2read,tmp] = uigetfile([Cfg.folder,'/*.xls*'],'Select organizer file to read');
if tmp==0, % cancel
    return
end
FileName = fullfile(Path2read,File2read);

[a b c] = xlsread(FileName);

% find keyword SEQUENCE
six = find(strcmp(c(:,1),'SEQUENCE'));

% loop though all found sequences
for k=1:length(six)
    
    % starting and finishing indexes
    n1 = six(k);
    if k<length(six), n2 = six(k+1)-1; else, n2 = length(c(:,1)); end
    
    % details
    S(k).filename = c(n1,2);
    S(k).time_start = datestr( cell2mat( c(n1+2,1) ) );
    S(k).env.temperature = cell2mat( c(n1+2,2) );
    S(k).env.humidity = cell2mat( c(n1+2,3) );
    S(k).env.wind_speed = cell2mat( c(n1+2,4) );
    S(k).env.noise_level_dba = cell2mat( c(n1+2,5) );
    
    S(k).calibration = cell2mat( c(n1+3,2:9) );
    S(k).distance = cell2mat( c(n1+4,2:9) );
    
    S(k).events.time = cell2mat( c(n1+6:n2,1) );
    S(k).events.duration = cell2mat( c(n1+6:n2,2) );
    S(k).events.type = c(n1+6:n2,3);
    
    for m=1:length(S(k).events.time)
        if isnan(S(k).events.time(m)) | isnan(S(k).events.duration(m)) 
            disp(sprintf('Error: sequence %d data organizer problem with event %d',k,m))
            return
        end
    end
end

N = length(S);

evtype=[];
for k=1:N
    seq_time(k,:) = S(k).time_start;
    temp(k)= S(k).env.temperature;
    humi(k)= S(k).env.humidity;
    wind(k)= S(k).env.wind_speed;
    noise(k)= S(k).env.noise_level_dba;
    calib(k,:)= S(k).calibration;
    evtype = [evtype;lower(S(k).events.type)];
end

figure(1)
subplot(221), plot(1:N,temp), set(gca,'xtick',1:N,'xticklabel',seq_time,'xticklabelrotation',45), ylabel('Temperature (\circF)'), axis tight, grid
subplot(222), plot(1:N,humi), set(gca,'xtick',1:N,'xticklabel',seq_time,'xticklabelrotation',45), ylabel('Humidity (%)'), axis tight, grid
subplot(223), plot(1:N,wind), set(gca,'xtick',1:N,'xticklabel',seq_time,'xticklabelrotation',45), ylabel('Wind speed (mph)'), axis tight, grid
subplot(224), plot(1:N,noise), set(gca,'xtick',1:N,'xticklabel',seq_time,'xticklabelrotation',45), ylabel('Noise (dBA)'), axis tight, grid

[ue, ia, ic] = unique(evtype);
for k=1:length(ue)
    nue(k) = length(find(ic==k));
end

set(handles.run_events,'string',ue);
set(handles.run_events,'value',1);

figure(2)
bar(1:length(ue), nue),set(gca,'xtick',1:length(ue),'xticklabel',ue,'xticklabelrotation',45), grid
ylabel('Events Occurence')

Cfg.unique_events.name = ue;
Cfg.unique_events.occurence = nue;

% figure(11), plot(1:N,calib(:,1:4),'+:'), set(gca,'ylim',[1e-7 6e-5]), legend('1','2','3','4')
% set(gca,'xtick',1:N,'xticklabel',seq_time,'xticklabelrotation',45), ylabel('Calibration factor (Pa/F8 a.u.)'),  grid

Cfg.org = S;
setappdata(CurFig,'Cfg',Cfg);



function event_number_Callback(hObject, eventdata, handles)
% hObject    handle to event_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of event_number as text
%        str2double(get(hObject,'String')) returns contents of event_number as a double
run_events_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function event_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to event_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton9_show_event.
function pushbutton9_show_event_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9_show_event (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'AudioViewerYDS');
Cfg = getappdata(CurFig,'Cfg');

FileNames = get(handles.listbox1_selected_data,'string'); 
EventValue = get(handles.run_events,'value');

events = cellstr(get(handles.run_events,'String')) ;
if isempty(events), return; end
selected_event = events{get(handles.run_events,'Value')};

%%
% check consistency
NumberOfFiles = length(FileNames);
NumberOfSequences = length(Cfg.org);

if NumberOfFiles~=NumberOfSequences,
    errordlg('Organizer: Number of files does not match the number of sequences');
    return
end
ok=0;
for k=1:NumberOfSequences
    if strcmp(char(Cfg.org(k).filename), char(FileNames{k}))
        disp(sprintf('Sequence %d matches file %s',k,char(FileNames{k})))
        ok=ok+1;
    else
        disp(sprintf('Sequence %d (%s) does not match file %s',k,char(Cfg.org(k).filename),char(FileNames{k})))
    end
end
if ok<NumberOfSequences
    errordlg('Organizer: There is a missmatch in organizer sequences definition and files');
    return
end

event_number = str2double(get(handles.event_number,'String'));
disp(sprintf('Processing event %s number %d',upper(selected_event),event_number))
n = event_number;

seqnum = get_sequence_num_for_event(handles, Cfg,selected_event,event_number);

chan_num = get(handles.show_channel,'value');
set(handles.show_all_channels,'value',0), drawnow;
set(handles.pressure_calibration_factor,'string',sprintf('%6.4e',Cfg.org(seqnum).calibration(chan_num)))

set(handles.show_start_time,'string',sprintf('%6.2f',Cfg.org(seqnum).events.time(event_number)))
set(handles.show_window_time,'string',sprintf('%6.2f',Cfg.org(seqnum).events.duration(event_number)))

filenum = get(handles.listbox1_selected_data,'Value');
if filenum ~= seqnum,
    set(handles.listbox1_selected_data,'Value',seqnum), drawnow;
    load_data_from_file(handles)
end

pushbutton3_show_Callback(hObject, eventdata, handles)

%%
% find sequence number for a given event
function seqnum = get_sequence_num_for_event(handles, Cfg,selected_event,event_number)
seqnum = [];
selev=[];
cnt = 1;

K = length(Cfg.org);

for k=1:K
    
    N = length(Cfg.org(k).events.type);

    for n=1:N
        if strcmp( lower(Cfg.org(k).events.type(n)),selected_event)
            selev(cnt)=k;
            cnt=cnt+1;
        end
    end
end

if isempty(selev), disp(sprintf('Error: unable to select event %s', upper(selected_event))), return, end
if length(selev)< event_number, disp(sprintf('Error: event number exceed occurence of event %s', upper(selected_event))), return, end

seqnum = selev(event_number);

str = sprintf('Sequence %d %s has selected event %d %s stating at %5.2f sec and lasting %5.2f seq',...
    seqnum, char(Cfg.org(seqnum).filename), event_number, upper(selected_event),...
    Cfg.org(seqnum).events.time(event_number),Cfg.org(seqnum).events.duration(event_number));

set(handles.events_text,'string',str)
disp(str)


% --- Executes on selection change in run_type.
function run_type_Callback(hObject, eventdata, handles)
% hObject    handle to run_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns run_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from run_type


% --- Executes during object creation, after setting all properties.
function run_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to run_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function down_sample_factor_Callback(hObject, eventdata, handles)
% hObject    handle to down_sample_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of down_sample_factor as text
%        str2double(get(hObject,'String')) returns contents of down_sample_factor as a double


% --- Executes during object creation, after setting all properties.
function down_sample_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to down_sample_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton10_filter.
function pushbutton10_filter_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'AudioViewerYDS');
Cfg = getappdata(CurFig,'Cfg');

N = 1; % this many spectral lines
for k=1:2*N,
    if mod(k,2)
        prompt{k} = sprintf('%02d start (kHz): ',k);
    else
        prompt{k} = sprintf('%02d stop (kHz): ',k);
    end
    defans{k} = 'NaN';
end
answer = inputdlg(prompt,'Filter Setup',1,defans);
filtfreq=[];
try
    for k=1:N
        n = 2*k-1;
        tmp1 = str2double(answer{n});
        tmp2 = str2double(answer{n+1});
        
        if ~isnan(tmp1) & ~isnan(tmp2)
            filtfreq(k,:) = [tmp1 tmp2]*1e3;
        end
    end
    if ~isempty(filtfreq)
        for k=1:length(filtfreq(:,1))
            n = 2*k-1;
            if filtfreq(n)>=filtfreq(n+1),
                disp('Error: filter window frequencies'), filtfreq=[]; break;
            end
            if k>1 
                if filtfreq(n)<filtfreq(n-1),
                    disp('Error: filter sequence frequencies'), filtfreq=[]; break;
                end
            end
            
        end
    end
    Cfg.filter_frequencies = filtfreq;
catch
    disp('Error: filter frequencies unrecognized input')
end
if ~isempty(Cfg.filter_frequencies)
    set(handles.filter_text,'string',sprintf('Filter set between %4.1f and %4.1f kHz',...
        Cfg.filter_frequencies(1)/1e3, Cfg.filter_frequencies(2)/1e3))
end
setappdata(CurFig,'Cfg',Cfg);
clc
pushbutton3_show_Callback(hObject, eventdata, handles);

function Sigout = FilterSingal(SigIn,Cfg,SamplingFrequency)
Sigout = SigIn;

if ~isempty(Cfg.filter_frequencies)
    fn1 = 2*Cfg.filter_frequencies(1)/SamplingFrequency;
    fn2 = 2*Cfg.filter_frequencies(2)/SamplingFrequency;
    if fn1<0 | fn1>=fn2 | fn1>=1, errordlg('Filter frequencies are not correct'), return, end
    if fn2>=1, [bb, aa] = butter( 4,fn1, 'high');
    else, [bb, aa] = butter( 4, [fn1 fn2] ); end
    for k=1:length(SigIn(1,:))
        Sigout(:,k) = filtfilt( bb, aa, SigIn(:,k) );
    end
end



function spl_db_Callback(hObject, eventdata, handles)
% hObject    handle to spl_db (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spl_db as text
%        str2double(get(hObject,'String')) returns contents of spl_db as a double


% --- Executes during object creation, after setting all properties.
function spl_db_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spl_db (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in test_signal_type.
function test_signal_type_Callback(hObject, eventdata, handles)
% hObject    handle to test_signal_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns test_signal_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from test_signal_type


% --- Executes during object creation, after setting all properties.
function test_signal_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to test_signal_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton12_generate_test_signal.
function pushbutton12_generate_test_signal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12_generate_test_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'AudioViewerYDS');
Cfg = getappdata(CurFig,'Cfg');
Data = getappdata(CurFig,'Data');

try
    Fs = Data.Fs;
catch
    Fs = 48000;
end


contents = cellstr(get(handles.test_signal_type,'String'));
sig_type = contents{get(handles.test_signal_type,'Value')}

switch sig_type
    case 'sin'
        answer = inputdlg({'Frequency (kHz):','Amplitude (Pa):','Duration (s):','Sampling frequency (kHz):','Bits per sample:'},...
            'Sin signal setup',1,{'1','1.0','1.0',sprintf('%6.4f',Fs),'24'});
%         answer = inputdlg({'Frequency (kHz):','Amplitude (Pa):','Duration (s):','Sampling frequency (kHz):',},...
%             'Sin signal setup',1,{'1','1.4142','10.0',sprintf('%6.4f',Fs)});
            f = 1000*str2num(answer{1});
            a = str2num(answer{2});
            t = str2num(answer{3});
            Fs = str2num(answer{4});
            Bps = str2num(answer{5});
            
            s = a*sin(2*pi*f*(0:(1/Fs):t));
            ScalingFactor = 2^Bps/2-1; 1/ScalingFactor
%             ScalingFactor = double(intmax('int32'))/100; 1/ScalingFactor
%             signal = signal./max(abs(signal(:)))*(1-(2^-(24-1)));
            s32 = int32(s*ScalingFactor); double([max(s32), min(s32)])/ScalingFactor
            audiowrite('audio_test_signal_sin.wav',s32,Fs,'BitsPerSample',Bps);
            
    otherwise,
        disp('TBD')
end



function mic_sens_Callback(hObject, eventdata, handles)
% hObject    handle to mic_sens (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mic_sens as text
%        str2double(get(hObject,'String')) returns contents of mic_sens as a double


% --- Executes during object creation, after setting all properties.
function mic_sens_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mic_sens (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sys_gain_Callback(hObject, eventdata, handles)
% hObject    handle to sys_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sys_gain as text
%        str2double(get(hObject,'String')) returns contents of sys_gain as a double


% --- Executes during object creation, after setting all properties.
function sys_gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sys_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton13_calibration_factor.
function pushbutton13_calibration_factor_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13_calibration_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurFig = getappdata(0,'AudioViewerYDS');
Cfg = getappdata(CurFig,'Cfg');
Data = getappdata(CurFig,'Data');

mic_sens = str2double(get(handles.mic_sens,'string'))/1000; % Pa/V
sys_gain = str2double(get(handles.sys_gain,'string'));

full_scale_dBV = -sys_gain + 15.17; % Zoom F8 station
full_scale_V = 10^(full_scale_dBV/20);

calibr_factor_Pa_per_sample = full_scale_V/2^23/mic_sens;
calibr_factor_Pa_per_V = full_scale_V/mic_sens;

disp(sprintf('Old calibration factor Pa/sample = %6.4e',calibr_factor_Pa_per_sample))
set(handles.pressure_calibration_factor,'string',sprintf('%6.4e',calibr_factor_Pa_per_V))


% --- Executes on selection change in spectrogram_type.
function spectrogram_type_Callback(hObject, eventdata, handles)
% hObject    handle to spectrogram_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns spectrogram_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from spectrogram_type
update_screen(handles)

% --- Executes during object creation, after setting all properties.
function spectrogram_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spectrogram_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
