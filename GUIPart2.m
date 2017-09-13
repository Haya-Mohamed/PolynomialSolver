function varargout = GUIPart2(varargin)
% GUIPART2 MATLAB code for GUIPart2.fig
%      GUIPART2, by itself, creates a new GUIPART2 or raises the existing
%      singleton*.
%
%      H = GUIPART2 returns the handle to a new GUIPART2 or the handle to
%      the existing singleton*.
%
%      GUIPART2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIPART2.M with the given input arguments.
%
%      GUIPART2('Property','Value',...) creates a new GUIPART2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUIPart2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUIPart2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUIPart2

% Last Modified by GUIDE v2.5 16-May-2017 23:22:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIPart2_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIPart2_OutputFcn, ...
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

function setGlobalval(x,y)
    global val
	global str
    val = x;
	str = y;
	
function r = getGlobalval
   global val
   r = val;	
	
function w = getGlobalStr
   global str
   w = str;
   

% --- Executes just before GUIPart2 is made visible.
function GUIPart2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUIPart2 (see VARARGIN)

% Choose default command line output for GUIPart2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUIPart2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUIPart2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
val = get(hObject, 'Value');
str = get(hObject, 'String');
setGlobalval(val,str);

switch str{val}
    case 'LU decomposition'
        set(handles.nPanel,'visible','on');
    case 'All Methods'
        set(handles.nPanel,'visible','on');
    otherwise
        set(handles.nPanel,'visible','off');
end

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iterBox_Callback(hObject, eventdata, handles)
% hObject    handle to iterBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iterBox as text
%        str2double(get(hObject,'String')) returns contents of iterBox as a double


% --- Executes during object creation, after setting all properties.
function iterBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iterBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function epsBox_Callback(hObject, eventdata, handles)
% hObject    handle to epsBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epsBox as text
%        str2double(get(hObject,'String')) returns contents of epsBox as a double


% --- Executes during object creation, after setting all properties.
function epsBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epsBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on InputTable and none of its controls.
function InputTable_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to InputTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when entered data in editable cell(s) in InputTable.
function InputTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to InputTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)



function nBox_Callback(hObject, eventdata, handles)
% hObject    handle to nBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nBox as text
%        str2double(get(hObject,'String')) returns contents of nBox as a double


% --- Executes during object creation, after setting all properties.
function nBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function inputArea_Callback(hObject, eventdata, handles)
% hObject    handle to inputArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputArea as text
%        str2double(get(hObject,'String')) returns contents of inputArea as a double


% --- Executes during object creation, after setting all properties.
function inputArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function inputPanel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputPanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SolveBtn.
function SolveBtn_Callback(hObject, eventdata, handles)
% hObject    handle to SolveBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%2*x + y + 3*z ==1, 2*x + 6*y + 8*z == 3, 6*x + 8*y+18*z==5
val = getGlobalval;
str = getGlobalStr;
eps = get(handles.epsBox,'String');
maxit = get(handles.iterBox,'String');
if(isempty(maxit))
    maxit = 50;
end
cla(handles.axes1);
equations = get(handles.inputArea,'String');
syms x y z
[A,B] = equationsToMatrix(equations, [x, y, z])
switch str{val}
    case 'Gaussian-elimination'
        [a,x,flag] = GaussElimination(A,B,0,handles)
        
        set(handles.outArea,'string',sprintf(strcat('GauessElim:','\n', char(x))));
    case 'LU decomposition'
        n = get(handles.nBox,'String');
        [x, er] = LUDecomposition(A, B, n, eps,0,handles)
        set(handles.outArea,'string',strcat('LU:',char(x)));
    case 'Gaussian-Jordan'
        matrix = cat(2, A, B);
        [ answers ] = guessJordon( matrix ,0,handles)
        set(handles.outArea,'string',strcat('GauessJordan:',char(answers)));
    case 'Gauss-Seidel'
        siz = size(B);
        x = zeros(siz);
        x = GaussSeidel(A,B,x,maxit,eps,0,handles);
        set(handles.outArea,'string',x);
        disp(x);
    case 'All Methods'
         n = get(handles.nBox,'String');
         siz = size(B);
         x = zeros(siz);
         matrix = cat(2, A, B);
% % % % % % % % % % % % % % % % % % % % % % %          
         disp('Gauss Elimination x')
         [a,x1,flag] = GaussElimination(A,B,1,handles)
         disp(x1);
% % % % % % % % % % % % % % % % % % % % % % %    
         disp('LU x')
         [x2, er] = LUDecomposition(A, B, n, eps,1,handles)
         disp(x2);
% % % % % % % % % % % % % % % % % % % % % % %    
         disp('Gauss Jordan x')
         [ answers ] = guessJordon( matrix,1,handles )
         disp(answers);
% % % % % % % % % % % % % % % % % % % % % % %  
         disp('Gauss Siedel x')
         z = GaussSeidel(A,B,x,maxit,eps,1,handles);
         disp(z);
         set(handles.outArea,'string',{sprintf(strcat('GauessElim:',char(x1),'\n','LU:','\n','GauessJordan:',char(answers),'\n','GauessSeidel:')),z});
end
% magneta
function x = GaussSeidel(A,B,x,maxit,ea,allMethods,handles)
flag = 1;
[na , ma ] = size (A);

if na ~= ma
   flag = 0;
    return
end
[nb , mb ] = size (B);
if nb ~= na || mb~=1
  flag = 0;
   return
end
 
[nx, mx] = size(x);
if nx ~= na || mx ~= 1
%     disp( 'ERROR: Check input')
%     return
end
iter = 0;
for r=1:maxit
    max = 0;
        for i=1: na
            s=0;
            for j=1: ma
               if j ~= i
                   s = s + (A(i,j)* x(j));
               end
            end
            t=(B(i)-s)/A(i,i);
             e = (abs(t-x(i)))/t;
             if e > max
                 max = e;
             end
            x(i)=t;
            error(iter+1, i) = e;
            if(allMethods == 1 && r+1<=i)
                axes(handles.axes1);
                line([r-1 r],[x(r) x(r+1)],'Color','m');
                hold on;
            end
        end
             if max < ea
                break;
             end 
        iter = r;
end
if iter > maxit
    flag = 0;
end
hold off;
% green  
function [a,x,flag] = GaussElimination(a,b,allMethods,handles)

flag = 1;
[na , ma ] = size(a);
if na ~= ma
    flag = 0;
    return
end
[nb , mb ] = size (b);
if nb ~= na || mb~=1
  flag =0;
   return
end

%forward elimination
   for k = 1:  size(a,1)-1
      c = a(k:size(a,1),k);
       [maxPivot,ind] = max(c(:))   
       a([k,ind+k-1],:) = a([ind+k-1, k],:);
       b([k,ind+k-1]) = b([ind+k-1, k]);
     for i = k+1:  size(a,1)
        factor = a(i,k) / a(k,k);
        for j = k : size(a,1)
             a(i,j) = a(i,j) - factor * a(k,j);
        end
        b(i) = b(i) - factor * b(k);
     end
   end
%back substitution
   x(size(a,1)) = b(size(a,1)) / a(size(a,1),size(a,1));
   cnt = 0;
   for i = size(a,1)-1: -1 : 1
      sum = 0;
      for j = i+1 : size(a,1)
         sum = sum + a(i,j) * x(j);
      end
      x(i) = (b(i) - sum) / a(i,i);
      if(allMethods == 1 && i > 1)
            axes(handles.axes1);
            line([cnt cnt+1],[x(i) x(i-1)],'Color','g');
            cnt = cnt + 1;
            hold on;
      end
   end
% red
function [ answers ] = guessJordon( matrix,allMethods,handles )
m=matrix;
flag=0;
for i=1:1:size(m,1)
    m=replace(m,i,maxJ(m,i));
    if(m(i,i)==0)
        flag = 1;
    end
    m(i,:)=m(i,:)/m(i,i);
    for r=1:1:size(m,1)
       if(r~=i)
       a=m(r,i);
       m(r,:)=m(r,:)-m(i,:)*a;
       if(allMethods == 1 && i > 1)
           axes(handles.axes1);
           line([r-1 r],[m(i) m(i+1)],'Color','r');
           hold on;
        end
       end
    end
end
disp(m);
if(flag==0)
answers = m(:,size(m,2));
else
    answers = NaN;
end

function [ Maxindex ] = maxJ( matrix , start)
maxI=start;
a=matrix(start,start);
for n=start:1:size(matrix,1)
 if(abs(matrix(n,start))>a)
     a=abs(matrix(n,start));
     maxI=n;
 end
end
Maxindex=maxI;

function [Matrix] = replace (matrix,x,y)
Martrixp=matrix;
a=Martrixp(x,:);
Martrixp(x,:)=Martrixp(y,:);
Martrixp(y,:)=a;
Matrix=Martrixp;
% blue
function [x, er] = LUDecomposition(a, b, n, tol,allMethods,handles)
    er = 0;
    n=size(b);
    s = zeros(n);
    o = zeros(n);
    for i = 1 : n 
        o(i) = i;
        s(i) = abs(a(i,1));
        for j = 2 : n
            if (abs(a(i,j)) > s(i))
                s(i) = abs(a(i,j));
            end
        end
    end
    for k = 1 : (n - 1) 
        p = k;
        big = abs(a(o(k), k) / s(o(k)));
        for i = k + 1 : n 
            dummy = abs(a(o(i), k) / s(o(i)));
            if (dummy > big) 
                big = dummy;
                p = i;
            end
        end
        dummy = o(p);
        o(p) = o(k);
        o(k) = dummy;
        if (abs(a(o(k),k) / s(o(k))) < tol) 
            er = -1;
            return;
        end
        for i = k + 1 : n 
            factor = a(o(i),k) / a(o(k),k);
            a(o(i),k) = factor;
            for j = k + 1 : n
                a(o(i), j) = a(o(i), j) - factor * a(o(k),j);
            end
        end
        for i = k + 1 : n
           s(i) = abs(a(i, 1));
            for j = k + 1 : n
                if (abs(a(i,j)) > s(i))
                    s(i) = abs(a(i,j));
                end
            end 
        end
    end
    if (abs(a(o(n), n) / s(o(n))) < tol)
        er = -1;
        return;
    end
    y = zeros([1 n]);
    y(o(1)) = b(o(1));
    for i = 2 : n 
        sum = b(o(i));
        for j = 1 : i - 1
            sum = sum - a(o(i), j) * y(o(j));
        end
        y(o(i)) = sum;
    end
    x = ones([1 n]);
    x(n) = y(o(n)) / a(o(n), n);
    for i = 1 : (n - 1) 
        sum = 0;
        for j = (n - i + 1) : n
            sum = sum + a(o(n - i), j) * x(j);
        end
        x(n - i) = (y(o(n - i)) - sum) / a(o( n - i), n - i);
        if(allMethods == 1 && i > 1)
           axes(handles.axes1);
           line([i-1 i],[x(n-i) x(n-i+1)],'Color','r');
           hold on;
        end
    end   


% --- Executes during object creation, after setting all properties.
function outArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
