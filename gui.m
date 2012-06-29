function varargout = gui(varargin)
% GUI MATLAB code for gui.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui

% Last Modified by GUIDE v2.5 29-Jun-2012 09:21:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_OutputFcn, ...
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


% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui (see VARARGIN)
contents = cellstr(get(handles.tree_name,'String'));
fn = contents{get(handles.tree_name,'Value')};
load(fn)
handles.Tree = Tree;
fprintf('Loaded %s \n', fn)
n_j = size(Tree,1);
js = cell(n_j, 1);
for i = 1:n_j,
    js{i} = num2str(i);
end
set(handles.j_level, 'String', js)
set(handles.j_level, 'Value', 1)
n_i = size(Tree{1,1}.ExtBasis, 2);
is = cell(n_i, 1);
for i = 1:n_i,
    is{i} = num2str(i);
end
set(handles.idx, 'String', is)
set(handles.idx, 'Value', 1)
do_plot(handles)

% Choose default command line output for gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in idx.
function idx_Callback(hObject, eventdata, handles)
% hObject    handle to idx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Tree = handles.Tree;
% j = get(handles.j_level, 'Value');
% i = get(hObject, 'Value');
do_plot(handles)

% Hints: contents = cellstr(get(hObject,'String')) returns idx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from idx


% --- Executes during object creation, after setting all properties.
function idx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to idx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in j_level.
function j_level_Callback(hObject, eventdata, handles)
% hObject    handle to j_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tree = handles.Tree;
j = get(hObject, 'Value');
s = 2 - get(handles.scaling, 'Value');
n_i = size(Tree{j,s}.ExtBasis, 2);
if n_i == 0,
    is = cell(1,1);
    is{1} = '---';
else
    is = cell(n_i, 1);
    for i = 1:n_i,
        is{i} = num2str(i);
    end
end
set(handles.idx, 'Value', 1)
set(handles.idx, 'String', is)

do_plot(handles)

% Hints: contents = cellstr(get(hObject,'String')) returns j_level contents as cell array
%        contents{get(hObject,'Value')} returns selected item from j_level


% --- Executes during object creation, after setting all properties.
function j_level_CreateFcn(hObject, eventdata, handles)
% hObject    handle to j_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tree_name.
function tree_name_Callback(hObject, eventdata, handles)
% hObject    handle to tree_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String'));
fn = contents{get(hObject,'Value')};
load(fn)
handles.Tree = Tree;
fprintf('Loaded %s \n', fn)
s = 2 - get(handles.scaling, 'Value');
n_j = size(Tree,s);
js = cell(n_j, 1);
for i = 1:n_j,
    js{i} = num2str(i);
end
set(handles.j_level, 'String', js)
set(handles.j_level, 'Value', 1)
n_i = size(Tree{1,s}.ExtBasis, 2);
is = cell(n_i, 1);
for i = 1:n_i,
    is{i} = num2str(i);
end
set(handles.idx, 'String', is)
set(handles.idx, 'Value', 1)

do_plot(handles)

guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns tree_name contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tree_name


% --- Executes during object creation, after setting all properties.
function tree_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tree_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

f = cell(4,1);
f{1} = 'circle_tree';
f{2} = 'Tree_grid_4';
f{3} = 'Tree_grid_8';
f{4} = 'Tree_grid_16';
set(hObject, 'String', f)
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in uipanel5.
function uipanel5_SelectionChangeFcn(hObject, eventdata, handles)
Tree = handles.Tree;

s = 2 - get(handles.scaling, 'Value');
j = get(handles.j_level, 'Value');
n_i = size(Tree{j,s}.ExtBasis, 2);
if n_i == 0,
    is = cell(1,1);
    is{1} = '---';
else
    is = cell(n_i, 1);
    for i = 1:n_i,
        is{i} = num2str(i);
    end
end
set(handles.idx, 'String', is)
set(handles.idx, 'Value', 1)
do_plot(handles)


function do_plot(handles)
Tree = handles.Tree;
s = 2 - get(handles.scaling, 'Value');
j = get(handles.j_level, 'Value');
if isempty(Tree{j,s}.ExtBasis)
    cla(handles.axes1)
    return
end
i = get(handles.idx, 'Value');
contents = cellstr(get(handles.tree_name,'String'));
name = lower(contents{get(handles.tree_name,'Value')});
x = Tree{j,s}.ExtBasis(:,i);
if strfind(name, 'grid')
    N = sqrt(length(x));
    imagesc(reshape(x,N,N))
else
    plot(x)
end
