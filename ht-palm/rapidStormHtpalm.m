function rapidStormHtpalm()
% RAPIDSTORMHTPALM - GUI for batch analysis of HTPALM
%TODO: delete the unneeded metadata files too (check which ones get used)
%TODO: add a warn dialog on the skip RS analysis option

%  Initialization tasks
%assign the default param values
param.dataFolder1 = '';
param.dataConfig1= '';
param.do3d=true;
param.calFile= 'G:\cameraSettings\astigmatic 3d calibration\150309 cal bead 0.1um -1.5to1.5 20nmstep.zCal.range-725To892.zCal.txt';
param.nDeleteSkip = 30;
param.doDeletePalmIm=true;
param.doDeleteOnly=false;
param.snr=50;
param.pixSzX=120;
param.pixSzY=120;
param.doTestRun=false;

%  Construct the components
dlg_title='HTPALM-RapidSTORM';
figPos = [360,500,280,340];
fh=figure('MenuBar','none','Name',dlg_title,'Toolbar','none','Position',figPos,'NumberTitle','off');

hStart = 35;
lAlign1 = 20;
lAlign2 = 30;
rAlign1 = 100;
%File path box
hCur = hStart;
uiCompHgt = 20;
ui.staticText_htpalmFile = uicontrol(fh,'Style','text',...
                'String','HTPALM config file:',...
                'Position',[lAlign1 figPos(4)-hCur 130 uiCompHgt]);
hCur = hCur+uiCompHgt +5;

ui.editBox_htpalmConfigFile = uicontrol(fh,'Style','edit',...
                'String','',...
                'Position',[lAlign2 figPos(4)-hCur 130 uiCompHgt]);
                
ui.pushButton_selectHtpalmConfigFile = uicontrol(fh,'Style','pushbutton','String','...',...
                'Position',[165 figPos(4)-hCur 30 uiCompHgt],...
                'Callback',@push_getPathHTPALM_conf);
hCur = hCur+uiCompHgt +5;
%TickBox Do3d Files
ui.checkBox_do3d = uicontrol(fh,'Style','checkbox',...
                'String','Do 3D fitting',...
                'Value',param.do3d,'Position',[lAlign1 figPos(4)-hCur 210 uiCompHgt]);
hCur = hCur+uiCompHgt +5;
%Cal path box
uiCompHgt = 20;
ui.staticText_calFile = uicontrol(fh,'Style','text',...
                'String','3D cal file/ 2D FWHM',...
                'Position',[lAlign1 figPos(4)-hCur 130 uiCompHgt]);
hCur = hCur+uiCompHgt +5;

ui.editBox_calFile = uicontrol(fh,'Style','edit',...
                'String',param.calFile,...
                'Position',[lAlign2 figPos(4)-hCur 130 uiCompHgt]);
ui.pushButton_selectCalFile = uicontrol(fh,'Style','pushbutton','String','...',...
                'Position',[165 figPos(4)-hCur 30 uiCompHgt],...
                'Callback',@push_getPathCalFile);
hCur = hCur+uiCompHgt +5;
%Param boxes
%SNR
ui.staticText_SNR = uicontrol(fh,'Style','text',...
                'String','SNR:',...
                'Position',[lAlign1 figPos(4)-hCur 70 uiCompHgt]); 
ui.editBox_SNR =  uicontrol(fh,'Style','edit',...
                'String',num2str(param.snr),...
                'Position',[rAlign1 figPos(4)-hCur 50 uiCompHgt]);
hCur = hCur+uiCompHgt +5;

%PixSz X
ui.staticText_pixX = uicontrol(fh,'Style','text',...
                'String','Pixel size X:',...
                'Position',[lAlign1 figPos(4)-hCur 70 uiCompHgt]); 
ui.editBox_pixX =  uicontrol(fh,'Style','edit',...
                'String',param.pixSzX,...
                'Position',[rAlign1 figPos(4)-hCur 50 uiCompHgt]);
hCur = hCur+uiCompHgt +5;
%PixSz Y
ui.staticText_pixY = uicontrol(fh,'Style','text',...
                'String','Pixel size Y:',...
                'Position',[lAlign1 figPos(4)-hCur 70 uiCompHgt]); 
ui.editBox_pixY =  uicontrol(fh,'Style','edit',...
                'String',param.pixSzY,...
                'Position',[rAlign1 figPos(4)-hCur 50 uiCompHgt]);
hCur = hCur+uiCompHgt +5;
%TickBox Delete Files
ui.checkBox_doDelete = uicontrol(fh,'Style','checkbox',...
                'String','Delete PALM images after processing',...
                'Value',param.doDeletePalmIm,'Position',[lAlign1 figPos(4)-hCur 210 uiCompHgt]);
hCur = hCur+uiCompHgt +5;
%NSkip
ui.staticText_nSkip = uicontrol(fh,'Style','text',...
                'String','Keep Nth file:','Value',param.doDeletePalmIm,...
                'Position',[lAlign1 figPos(4)-hCur 70 uiCompHgt]); 
ui.editBox_nSkip=  uicontrol(fh,'Style','edit',...
                'String',num2str(param.nDeleteSkip),...
                'Position',[rAlign1 figPos(4)-hCur 50 uiCompHgt]);
%Delete only tickbox
ui.checkBox_deleteOnly= uicontrol(fh,'Style','checkbox',...
                'String','Skip RapidStorm (delete only)',...
                'Value',param.doDeleteOnly,'Position',[lAlign1 figPos(4)-hCur 210 uiCompHgt]);
hCur = hCur+uiCompHgt +5;
%Do test run tickbox
ui.checkBox_doTestRun= uicontrol(fh,'Style','checkbox',...
                'String','Test Run (first movie only)',...
                'Value',param.doTestRun,'Position',[lAlign1 figPos(4)-hCur 210 uiCompHgt]);
hCur = hCur+uiCompHgt +5;
% Ok, cancel buttons
hCur = hCur+uiCompHgt +5;
uiCompHgt =30;
ui.pushbutton_ok = uicontrol(fh,'Style','pushbutton','String','Ok',...
                'Position',[35 figPos(4)-hCur 50 uiCompHgt],...
                'Callback', @pushOk);
ui.pushbutton_cancel = uicontrol(fh,'Style','pushbutton','String','Cancel',...
                'Position',[150 figPos(4)-hCur 50 uiCompHgt],...
                'Callback', @pushCancel);

%  Callbacks
  function pushOk(hObject,eventdata)
    if ~isempty(get(ui.editBox_htpalmConfigFile,'String'))&&~isempty(get(ui.editBox_calFile,'Value'));
      [pathstr, name, ext]= fileparts(get(ui.editBox_htpalmConfigFile,'String'));
      param.dataFolder1 = pathstr;
      param.dataConfig1 = [name,ext];

      param.do3d = get(ui.checkBox_do3d,'Value');
      if param.do3d
        param.calFile = get(ui.editBox_calFile,'String');
        param.fwhm=NaN;
      else
        param.calFile= '';
        param.fwhm = str2num(get(ui.editBox_calFile,'String'));
      end
      param.nDeleteSkip = str2num(get(ui.editBox_nSkip,'String'));
      param.doDeletePalmIm=get(ui.checkBox_doDelete,'Value');
      param.doDeleteOnly=get(ui.checkBox_deleteOnly,'Value');
      param.doTestRun = get(ui.checkBox_doTestRun,'Value');
      param.snr= str2num(get(ui.editBox_SNR,'String'));
      param.pixSzX= str2num(get(ui.editBox_pixX,'String'));
      param.pixSzY= str2num(get(ui.editBox_pixY,'String'));

      close(fh);
      runHtpalmRS(param);
    else
      errordlg('Config file or calibration file not chosen!');
    end
  end

  function pushCancel(hObject,eventdata)
    close(fh);
  end

  function push_getPathHTPALM_conf(hObject,eventdata)
    [FileName,PathName] = uigetfile('*config.xml', 'Select the HTPALM config file');
    set(ui.editBox_htpalmConfigFile,'String',fullfile(PathName,FileName));
  end

  function push_getPathCalFile(hObject,eventdata)
    [FileName,PathName] = uigetfile('*.zCal.txt', 'Select the Z-calibration file', param.calFile);
    set(ui.editBox_calFile,'String',fullfile(PathName,FileName));
  end
end
%------------------------
function runHtpalmRS(param)
param
doRun = true;
tformSaveName='';
mTrackParamFile='';
zScaleFactor=[];
xyWobbleFile='';
finalRoiBorder='';
srPixSize = [param.pixSzX, param.pixSzY];
testRun = false;

htpalmArg ={'RapidStormOnly'};
if param.doDeletePalmIm
  htpalmArg ={htpalmArg{:},'DeletePalm','NDeleteSkip',param.nDeleteSkip};
end

if param.doDeleteOnly
  htpalmArg ={htpalmArg{:},'SkipRapidStorm'};
end

if param.doTestRun 
  htpalmArg ={htpalmArg{:},'TestRun'};
end

if ~param.do3d
  htpalmArg = {htpalmArg{:},'2DFit',param.fwhm};
end
htpalmArg
if param.doDeletePalmIm && param.doDeleteOnly
  choice = questdlg('Pressing OK will delete PALM images WITHOUT analysing them!','Warning','Ok','Cancel','Cancel');
  if ~strcmp(choice,'Ok')
    doRun=false;
  end
end

if doRun
  analyseHtpalm(param.dataFolder1,param.dataConfig1,tformSaveName,srPixSize, param.snr, param.calFile,mTrackParamFile,zScaleFactor,xyWobbleFile,finalRoiBorder,htpalmArg{:});
end

end

