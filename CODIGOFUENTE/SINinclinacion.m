function [FSD,FSV,FSD1,FSV1,FSD2,FSV2,FSD3,FSV3]= SINinclinacion(Base,Ha,AlTalo,AncT,Co,concreto,So,Cant,d,Fp,Mp,hcopiaa,k,phii,phi,c,nf,n,R,h,g,R1,R2,BaseP,aa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLANO
                 if ismcc || isdeployed
                    % Make sure DOM is compilable
                    makeDOMCompilable();
                 end
                import mlreportgen.dom.*
                import mlreportgen.report.*
            
                assignin('base','PHI',phi);
                %condicional de cambio de estrato ampliar matriz de alturas
                 hcopia=h;
                if nf > 0
                    h(n+1)=h(n);
                end
            
            
                %hallar ka o ko rankie plano
                if k==2
                    ka=zeros(n,1);
                    for j=1:n
                        ka(j)=(1-sind(phii(j)))/(1+sind(phii(j)));
                    end
                    assignin('base','Ka',ka);
                end
                if k==3 %cuando es Ko reposo
                    ka=zeros(n,1);
                    for j=1:n
                        ka(j)=1-sind(phii(j));
                    end
                     assignin('base','Ko',ka);
                end
                %se agrega pasivo 
                if k==4
                    ka=zeros(n,1);
                    for j=1:n
                        ka(j)=(1+sind(phii(j)))/(1-sind(phii(j)));
                    end
                    assignin('base','Ka',ka);
                end

                ka=round(ka,3);
            
            
            
            
            
                %Condicional si esta humedo o no
                if nf>0
                    if nf~=1
                        %antes del nivel freático
                        for kk=1:nf-1
                            sv(kk)=h(kk)*g(kk);
                        end
                    end
                    %después del nivel freático
                    for kk2=nf:n
            
                        if R1==1
                        sv(kk2)=h(kk2)*(g(kk2)-9.81);
                        end
                        if R2==1
                            sv(kk2)=h(kk2)*(g(kk2)-62.4);
                        end
                    end
                end
                if nf==0
                    %ESTA SECO
                    for kk=1:n
                        sv(kk)=h(kk)*g(kk);
                    end
            
                end
                assignin('base','Sv',sv);
                %Carga
                for kk4=1:n
                    Soo(kk4)=So*ka(kk4);
                end
                assignin('base','SOO',Soo);
            
                for kk5=1:n
                    PesoCarga(kk5)=Soo(kk5)*h(kk5);
                end
                assignin('base','Pcarga',PesoCarga);
                Pcarga=sum(PesoCarga);
                assignin('base','PCARG',Pcarga);
            
                %esfuerzo total
                St=0;
            
                for kk3=2:n+1
                    St(kk3)=St(kk3-1)+sv(kk3-1);
                end
                 assignin('base','St',St);
            
            
                %%%%%%%%%%%%%%%%%%%%%%%%%%%esfuerzo horizontal
                %%%%%%%%%%%esfuerzo horizontal
            
            
            
                % Caso donde phis difernetes para estrato 1 y 2
                if n==1 || phi(1)~=phi(2)
                    for ii=1:n
                        sa(ii)=St(ii)*ka(ii);
                        sb(ii)=St(ii+1)*ka(ii);
                        assignin('base','Sb',sb);
                        assignin('base','Sa',sa);
                    end
                    %reorganiza
                    w=1;
                    for i=1:n
                        SS(w)=sa(i);
                        w=w+1;
                        SS(w)=sb(i);
                        w=w+1;
                        assignin('base','SS',SS);
                    end
            
                    %Fuerzas actuantes
            
            
                    if n==1
                        F(1)=h(1)*SS(2)/2;
                    end
            
                    %%%oswa
                    if n>1
                        %Reorganizar tamaño de h
                        sizeh=size(hcopia);
                        assignin('base','hcopia',hcopia);
            
                        sizeh=sizeh(1);
                        assignin('base','sizeh',sizeh);
            
                        %CAMBIE PARA 2 ESTRATOS PARA QUE FUNCIONARA
                        if n==2
                            h0=zeros(1,sizeh+1);
                        end
                        if n>2
                             h0=zeros(1,sizeh+2);
                        end
            
            
                        sizeh0=size(h0);
                        sizeh0=sizeh0(2);
                        h0(1)=h(1);
                        pos=2;
            
                        for i=2:sizeh
                            h0(pos)= h(i);
                            pos= pos+1;
                            h0(pos) = h(i);
                            pos= pos+1;
            
                        end
                         assignin('base','H0',h0);
            
            
                        F=zeros(1,sizeh0);
                        F(1)=h(1)*SS(2)/2;
            
                        for i=2:sizeh0
                            residuo=rem(i,2);
                            if residuo==0
                                F(i)=h0(i)*SS(i+1);
                            else
                                F(i)=h0(i)*(SS(i+1)-SS(i))/2;
                            end
            
            
                        end
            
            
                        assignin('base','Factua',F);
                    end
                end
            
            
                %%%%%%%%%%%%%%%%%%%% phis iguales en estratos 1 y 2
                if n>=2
            
                    if ka(1)== ka(2)%si phi son iguales cambia el tamaño de las matrices de fuerzas
            
                        for ii=1:n
                            sa(ii)=St(ii)*ka(ii);
                        end
                        if n==2
                            sb(1)=St(2)*ka(1);
                            sb(2)=St(3)*ka(2);
                            %reorganiza
                            w=1;
                            for i=1:n
                                SS(1)=sa(1);
                                w=w+1;
                                SS(w)=sb(i);
                            end
                            F(1)=SS(3)*(h(1)+h(2))/2;
                            assignin('base','Factua',F);
                        end
            
                        if n==3
                            sb(1)=St(3)*ka(1);
                            sb(2)=St(3)*ka(3);
                            sb(3)=St(4)*ka(3);
                            %reorganiza
                            w=1;
                            for i=1:n
                                SS(1)=sa(1);
                                w=w+1;
                                SS(w)=sb(i);
            
                            end
            
                            F(1)=SS(2)*(h(1)+h(2))/2;
                            F(2)=SS(3)*(h(3));
                            F(3)=(SS(4)-SS(3))*(h(3))/2;
                            assignin('base','Factua',F);
                        end
                        if n==4
                            sb(1)=St(3)*ka(1);
                            sb(2)=St(3)*ka(3);
                            sb(3)=St(4)*ka(3);
                            sb(4)=St(5)*ka(4);
                            %reorganiza
                            w=1;
                            for i=1:n
                                SS(1)=sa(1);
                                w=w+1;
                                SS(w)=sb(i);
            
                            end
                            F(1)=SS(2)*(h(1)+h(2))/2;
                            F(2)=SS(3)*(h(3));
                            F(3)=(SS(4)-SS(3))*(h(3))/2;
                            F(4)=SS(5)*(h(4));
                            F(5)=(SS(5)-SS(4))*h(4)/2;
                            assignin('base','Factua',F);
                        end
            
            
                    end
                end
            
            
            
                 assignin('base','Sa',sa);
                 assignin('base','SB',sb);
                 assignin('base','SS',SS);
            
                % Peso agua si hay nivel freatico
            
                if nf>0
                    %Nivel freatico en 1 y solo hay dos estratos
                    if nf==1 && n==2
                        if R1==1
                            pesoagua=9.81*(h(1)+h(2))^2/2;
                            assignin('base','PesoAgua',pesoagua);
                        end
                        if R2==1
                            pesoagua=62.4*(h(1)+h(2))^2/2;
                            assignin('base','PesoAgua',pesoagua);
                        end
                    end
                    %Nivel freatico mayor que uno y hay dos estratos
                    if nf>1 && n==2
                        if R1==1
                            pesoagua=9.81*(h(nf))^2/2;
                            assignin('base','PesoAgua',pesoagua);
                        end
                        if R2==1
                            pesoagua=62.4*(h(nf))^2/2;
                            assignin('base','PesoAgua',pesoagua);
                        end
                    end
                    %nivel freatico igual a uno y hay mas de dos estratos
                    if n>2 && nf==1
                        if R1==1
                            pesoagua=9.81*(sum(h))^2/2;
                            assignin('base','PesoAgua',pesoagua);
                        end
                        if R2==1
                            pesoagua=62.4*(sum(h))^2/2;
                            assignin('base','PesoAgua',pesoagua);
                        end
                    end
            
                    %Nivel freatico si es mayor igual que dos y hay mas de dos estratos
                    if nf>=2 && n>2
                        hf=0;
                        for i=nf:n
                            hf=hf+h(i);
                        end
                        if R1==1
                            pesoagua=9.81*(hf)^2/2;
                        end
                        if R2==1
                            pesoagua=62.4*(hf)^2/2;
                        end
                    end
            
                else
                    pesoagua=0;
                end
            
            
                PesoAct=sum(F)+pesoagua;
                assignin('base','PeACT',PesoAct);
            
                %%%%%%%%%Resistente
                %Volumen
                if Cant==1
                    Vol(1)=Co*(Ha-AlTalo)*1;
                    Vol(2)=(BaseP-Co)*(Ha-AlTalo)/2;%Parte triangular
                    Vol(3)=Base*AlTalo*1;
                    hactuante=h;
                    hactuante(n)=hactuante(n)-AlTalo;
                    for i=1:n
                    Vol(i+3)=hactuante(i)*AncT;
                    end
                    assignin('base','VOLUM',Vol);
                    %Peso
                    Pes(1)=Vol(1)*concreto;%Kn
                    Pes(2)=Vol(2)*concreto;%Kn
                    Pes(3)=Vol(3)*concreto;%Kn
                    for i=1:n
                    Pes(i+3)=Vol(i+3)*g(i);
                    end
                end
            
                if Cant==0
                    Vol(1)=Co*(Ha-AlTalo)*1;
                    Vol(2)=Base*AlTalo*1;
                    hactuante=h;
                    hactuante(n)=hactuante(n)-AlTalo;
                    for i=1:n
                    Vol(i+2)=hactuante(i)*AncT;
                    end
                    assignin('base','VOLUM',Vol);
                    %Peso
                    Pes(1)=Vol(1)*concreto;%Kn
                    Pes(2)=Vol(2)*concreto;%Kn
                    for i=1:n
                    Pes(i+2)=Vol(i+2)*g(i);
                    end
                end
            
            
                assignin('base','PESO',Pes);
            
                %PESO SIN CARGA
                PesoRes=sum(Pes);
                assignin('base','PesoSin',PesoRes);
            
               %Brazos Resistente
            
               anctalon=Base-AncT-Co;
            
               if Cant==1
                   B(1)=anctalon+(Co/2);
                   B(2)=anctalon+(Co*2/3);
                   B(3)=Base/2;
                   for i=1:n
                   B(i+3)=(AncT/2)+Co+anctalon;
                   end
                   assignin('base','Brazos',B);
            
                   %Momento resistente
                   M(1)=B(1)*Pes(1);
                   M(2)=B(2)*Pes(2);
                   M(3)=B(3)*Pes(3);
                   for i=1:n
                   M(i+3)=B(i+3)*Pes(i+3);%desde 4 es suelo
                   end
               end
               if Cant==0
                   B(1)=anctalon+(Co/2);
                   B(2)=Base/2;
                   for i=1:n
                   B(i+2)=(AncT/2)+Co+anctalon;
                   end
                   assignin('base','Brazos',B);
                   %Momento resistente
                   M(1)=B(1)*Pes(1);
                   M(2)=B(2)*Pes(2);
            
                   for i=1:n
                       M(i+2)=B(i+2)*Pes(i+2);%desde 3 es suelo
                   end
            
               end
            
            
               Mresistente=sum(M);
            
               %%%%%%%%%%%%%%%Brazo Actuante
               br(1)=h(1)*2/3;
               j=2;
               Lr=h(1);
            
                %Mayor a uno con phi diferente
               if n>1 && ka(1)~=ka(2)
                   %brazos de arriba hacia abajo
                   for i=2:n
                       br(j)=Lr+h(i)/2;
                       br(j+1)=Lr+h(i)*2/3;
                       j=j+2;
                       Lr=Lr+h(i);
                   end
               end
            
               %Si solo hay un estrato
               if n==1
                   br(1)=h(1)*2/3;
               end
               %mas de 2 estratos con k1 y k2 iguales salen 3 brazos si son 3 estratos
               if n>2 && ka(1)== ka(2)
                   hka=h(1)+h(2);
                   Lr=hka;
                   br(1)=hka*2/3;
                   if n==3
                       br(2)=Lr+h(3)/2;
                       br(3)=Lr+h(3)*2/3;
                   end
                   if n==4
                       br(2)=Lr+h(3)/2;
                       br(3)=Lr+h(3)*2/3;
                       Lr=Lr+h(3);
                       br(4)=Lr+h(4)/2;
                       br(5)=Lr+h(4)*2/3;
            
                   end
            
               end
            
            
            
               assignin('base','br',br);
               %brazos en el punto de abajo
               brvolteados=sum(hcopiaa)-br;
               assignin('base','Brazosr',brvolteados);
               %%%brazo agua
               if nf==0
                   brazoagua=0;
               else
                   brazoagua=h(nf)/3;

               end
               
               magua=brazoagua*pesoagua;
               
               %Momento actuante
               Maa=brvolteados.*F;
               assignin('base','Maa',Maa);
               Mactuante=sum(Maa)+magua;
               assignin('base','Mactuante',Mactuante);
            
               %%%%%Brazos carga
               bc(1)=h(1)*1/2;
               k=2;
               Lrr=h(1);
               %brazos de arriba hacia abajo
               for i=2:n
                   bc(k)=Lrr+h(i)/2;
                   k=k+1;
                   Lrr=Lrr+h(i);
               end
               assignin('base','BrCarga',bc);
               bcvolteados=sum(hcopiaa)-bc;
               assignin('base','bc',bcvolteados);
               %Momentos Carga
               Mc=PesoCarga.*bcvolteados;
               Mcarga=sum(Mc);
            
                %%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%R1=1
                %Caso con suelo con posterior
                if R==1
                    FSD=(PesoRes*tand(d))/(PesoAct-Fp);
                    FSV=Mresistente/(Mactuante);
            
%                     if So==0
%                          set(handles.FSD,'String',round(FSD,2));
%                 %        set(handles.FSDc,'String',FSDc);
%                          set(handles.FSV,'String',round(FSV,2));
%             
%                     end
                    
                    %Caso con carga despues del talon
                    FSD2=(PesoRes)*tand(d)/((PesoAct-Fp)+Pcarga);
                    FSV2=Mresistente/((Mactuante)+Mcarga);
                    assignin('base','FSD2',FSD2);
                    assignin('base','FSV2',FSV2);
%                     if aa==2
%                         set(handles.FSD,'String',round(FSD1,2));
%                         set(handles.FSV,'String',round(FSV1,2));
%                     end
                    %Caso con carga antes del talon
                    FSD1=(PesoRes+Pcarga)*tand(d)/(PesoAct-Fp);
                    FSV1=(Mresistente+Mcarga)/(Mactuante);
                    assignin('base','FSD1',FSD1);
                    assignin('base','FSV1',FSV1);
                    %Caso con carga antes y despues  del talon
                    FSD3=(PesoRes+Pcarga)*tand(d)/((PesoAct-Fp)+Pcarga);
                    FSV3=(Mresistente+Mcarga)/((Mactuante)+Mcarga);
                     assignin('base','FSD3',FSD3);
                     assignin('base','FSV3',FSV3);
                    
                end
            
                %%%%%%%%%%%%Caso sin suelo posterior posterior%%%%%%R=0
                if R==0
            
                    FSD=PesoRes*tand(d)/PesoAct;
                    %FSDc=(PesoRes*tand(d)+(c1*Base))/PesoAct;
                    FSV=Mresistente/Mactuante;
                    if So==0
            
                    
                    end
            
            
                        %Caso con carga despues del talon
                        FSD2=PesoRes*tand(d)/(PesoAct+Pcarga);
                        FSV2=Mresistente/(Mactuante+Mcarga);
                        assignin('base','FSD2',FSD2);
                        assignin('base','FSV2',FSV2);
            
                        %Caso con carga antes del talon
                        FSD1=(PesoRes+Pcarga)*tand(d)/PesoAct;
                        FSV1=(Mresistente+Mcarga)/Mactuante;
                        assignin('base','FSD1',FSD1);
                        assignin('base','FSV1',FSV1);
                        %Caso con carga antes y despues  del talon
                        FSD3=(PesoRes+Pcarga)*tand(d)/(PesoAct+Pcarga);
                        FSV3=(Mresistente+Mcarga)/(Mactuante+Mcarga);
                         assignin('base','FSD3',FSD3);
                         assignin('base','FSV3',FSV3);
            
                %         set(handles.cargaantes,'String',FSD2);
                %         set(handles.cargadespues,'String',FSD1);
                %         set(handles.cargaambos,'String',FSD3);
                %
                %         set(handles.maantes,'String',FSV2);
                %         set(handles.madespues,'String',FSV1);
                %         set(handles.maambos,'String',FSV3);
            
                    assignin('base','FSDD',FSD);
                    assignin('base','FSvv',FSV);
                end
            %%%%%%%%%%%%%%%GENERAR REPORTE
            if n==2
                reporte=questdlg('Desea generar reporte?', ...
	                        'Generar reporte', ...
	                        'SI',' NO  ','  CANCELAR  ','');
                        % Handle response
                        switch reporte
                            case 'SI' 
                                                             
           

                                % Add content to container (required)
                                % Types of content added here: title 
                                % page and table of contents reporters
                                doctype = "pdf";                  
                                rpt = Report("REPORTE MUROS DE RETENCIÓN",doctype);
                                open(rpt);
                                
                    
                                layoutObj = rpt.Layout;
                                layoutObj.PageSize = PageSize("10in","11in","portrait");
                    
                                %Create Report container
                                rpt.Layout.Landscape = 1;
                                rpt.Layout.PageNumberFormat = "i";
                                titlepg = TitlePage;
                    
                                 %%%%%%%%%%%%
                                titlepg.Title = 'ANÁLISIS GEOTÉCNICO DE MUROS DE RETENCIÓN';
                                titlepg.Author = 'Juan Solarte';
                                tplayoutObj = titlepg.Layout;
                                tplayoutObj.PageSize = PageSize("10in","12in","portrait");
                                imDir = which("m.png");
                                if exist(imDir,"file") == 2
                                    ImObj = Image(imDir);
                                    ImObj.Height = "5in";
                                    ImObj.Width = "6in";
                                    titlepg.Image = ImObj;
                                end
                                add(rpt,titlepg);
                                add(rpt,TableOfContents);
                        
                    %%%%%%%%%%%%%%%5
                                % Add content to report sections (optional)
                                % Text and formal image added to chapter
                    
                    
                                chap = Chapter('Valores de entrada');
                                %TABLA MURO Triangular
                                if Cant==0


                                    s1=Section("CARACTERISTICAS MURO(Triangular)");
                                    seccion1=["ALTURA MURO","BASE","CORONA ","ANCHO TALÓN","ALTURA TALÓN","PESO ESPECIFICO (concreto)"]';
                                    volumen1=[Ha,Base,Co,AncT,AlTalo,concreto]';
                                   
                                    
                                    Tabla1= FormalTable(["PROPIEDAD","VALOR"],[seccion1,volumen1]);
                                    Tabla1.Header.Style{end+1} = BackgroundColor("silver");
                                    Tabla1.Width = "500pt";
                                    Tabla1.Border="solid";
                                    Tabla1.ColSep="solid";
                                    Tabla1.RowSep="solid";
                                    Tabla1.HAlign="left";
                                    Tabla1.TableEntriesHAlign="center";

                                end


                                 %Tabla muro canteliver

                                if Cant==1
                                    %Base,Ha,AlTalo,AncT,Co,concreto,So,Cant,d,Fp,Mp,hcopiaa,k,phii,phi,c,nf,n,R,h,g,R1,R2,BaseP)

                                    s1=Section("CARACTERISTICAS MURO(Canteliver)");
                                    seccion1=["ALTURA MURO","BASE","CORONA ","ANCHO TALÓN","ALTURA TALÓN","BASE PANTALLA","PESO ESPECIFICO (concreto)"]';
                                    volumen1=[Ha,Base,Co,AncT,AlTalo,BaseP,concreto]';
                                   
                                    
                                    Tabla1= FormalTable(["PROPIEDAD","VALOR"],[seccion1,volumen1]);
                                    Tabla1.Header.Style{end+1} = BackgroundColor("silver");
                                    Tabla1.Width = "500pt";
                                    Tabla1.Border="solid";
                                    Tabla1.ColSep="solid";
                                    Tabla1.RowSep="solid";
                                    Tabla1.HAlign="left";
                                    Tabla1.TableEntriesHAlign="center";

                                    
                                
                                end

                               
                                
                                %TABLA SUELO RELLENO
                                 ss1=Section("CARACTERISTICAS SUELO DE RELLENO");
                                  %Base,Ha,AlTalo,AncT,Co,concreto,So,Cant,d,Fp,Mp,hcopiaa,k,phii,phi,c,nf,n,R,h,g,R1,R2,BaseP)
                                
                                seccions1=["No. ESTRATOS","ESTRATO DEL NIVEL FREATICO" ]';                           
                                                                                                                                   
                                volumens1=[n,nf]';                              
                                  
                                
                                Tablas1= FormalTable(["PROPIEDAD","VALOR"],[seccions1,volumens1]);
                                Tablas1.Header.Style{end+1} = BackgroundColor("silver");
                                Tablas1.Width = "300pt";
                                Tablas1.Border="solid";
                                Tablas1.ColSep="solid";
                                Tablas1.RowSep="solid";
                                Tablas1.HAlign="left";
                                Tablas1.TableEntriesHAlign="center";


                                %%%
                                seccions11=[1,2]';                         
                                                                                                                                   
                                volumens11=[hcopia(1),hcopia(2)]'; 
                                p=[g(1),g(2)]'; 
                                pi=[phi(1),phi(2)]'; 
                                if So==0
                                    Tablas11= FormalTable(["ESTRATO","ALTURA ESTRATO","PESO ESPECIFICO ESTRATO","ANGULO PHI ESTRATO"],[seccions11,volumens11,p,pi]);
                                    Tablas11.Header.Style{end+1} = BackgroundColor("silver");
                                    Tablas11.Width = "500pt";
                                    Tablas11.Border="solid";
                                    Tablas11.ColSep="solid";
                                    Tablas11.RowSep="solid";
                                    Tablas11.HAlign="left";
                                    Tablas11.TableEntriesHAlign="center";
                                    Tablas11.Style= [Tablas11.Style {NumberFormat("%.2f")}];
                                end
                                
                                %si hay una carga
                                if So>0
                                    if aa==2
                                        sor=[So,""]';
                                        Tablas11= FormalTable(["ESTRATO","ALTURA ESTRATO","PESO ESPECIFICO ESTRATO","ANGULO PHI ESTRATO","SOBRECARGA(antes del talón)"],[seccions11,volumens11,p,pi,sor]);
                                        Tablas11.Header.Style{end+1} = BackgroundColor("silver");
                                        Tablas11.Width = "500pt";
                                        Tablas11.Border="solid";
                                        Tablas11.ColSep="solid";
                                        Tablas11.RowSep="solid";
                                        Tablas11.HAlign="left";
                                        Tablas11.TableEntriesHAlign="center";
                                        Tablas11.Style= [Tablas11.Style {NumberFormat("%.2f")}];
                                    end

                                    if aa==3
                                        sor=[So,""]';
                                        Tablas11= FormalTable(["ESTRATO","ALTURA ESTRATO","PESO ESPECIFICO ESTRATO","ANGULO PHI ESTRATO","SOBRECARGA(despues del talón)"],[seccions11,volumens11,p,pi,sor]);
                                        Tablas11.Header.Style{end+1} = BackgroundColor("silver");
                                        Tablas11.Width = "500pt";
                                        Tablas11.Border="solid";
                                        Tablas11.ColSep="solid";
                                        Tablas11.RowSep="solid";
                                        Tablas11.HAlign="left";
                                        Tablas11.TableEntriesHAlign="center";
                                        Tablas11.Style= [Tablas11.Style {NumberFormat("%.2f")}];
                                    end
                                    if aa==4
                                        sor=[So,""]';
                                        Tablas11= FormalTable(["ESTRATO","ALTURA ESTRATO","PESO ESPECIFICO ESTRATO","ANGULO PHI ESTRATO","SOBRECARGA(antes y despues del talón)"],[seccions11,volumens11,p,pi,sor]);
                                        Tablas11.Header.Style{end+1} = BackgroundColor("silver");
                                        Tablas11.Width = "500pt";
                                        Tablas11.Border="solid";
                                        Tablas11.ColSep="solid";
                                        Tablas11.RowSep="solid";
                                        Tablas11.HAlign="left";
                                        Tablas11.TableEntriesHAlign="center";
                                        Tablas11.Style= [Tablas11.Style {NumberFormat("%.2f")}];
                                    end
                                    

                                end

                               

                                

                                %%%
                    
                                %TABLA SUELO DE CIMENTACION
                                ssr1=Section("CARACTERISTICAS SUELO DE CIMENTACIÓN");
                                seccionsr1=["ANGULO DE FRICCION Base (ϕ2)","PESO ESPECIFICO (cimentación)","COHESIÓN"]';
                                volumensr1=[phii(end),g(end),c]';
                               %phib,concreto,espesor,phi,PesoE,Cohesion,n,k,espe,pe,ph,R,TS,R1,R2,kh,kv,So
                                
                                
                                Tablasr1= FormalTable(["PROPIEDAD","VALOR"],[seccionsr1,volumensr1]);
                                Tablasr1.Header.Style{end+1} = BackgroundColor("silver");
                                Tablasr1.Width = "500pt";
                                Tablasr1.Border="solid";
                                Tablasr1.ColSep="solid";
                                Tablasr1.RowSep="solid";
                                Tablasr1.HAlign="left";
                                Tablasr1.TableEntriesHAlign="center";
                                Tablasr1.Style= [Tablasr1.Style {NumberFormat("%.2f")}];
                                
                                add(s1,Tabla1)
                               
                                add(ss1,Tablas1)
                                add(ss1,Tablas11)
                                add(ssr1,Tablasr1)
                                add(chap,s1)
                                add(chap,ss1)
                                add(chap,ssr1)
                                add(rpt,chap)
                    
                                 %%%%%%%Capitulo2
                                chap2 = Chapter('Factores Resistentes');
                                s=Section("Partes resistentes");

                                if Cant==0
                                    image2 = mlreportgen.report.FormalImage();
                                    image2.Image = which('reportetriaes.png');
                                    image2.Caption = 'Partes Resistente(Muro Triangular)';
                                    image2.Height = '3in';
                                    image2.Width = '4in';   

                                    seccion=["A","B","C","D"," "]';
                                    volumen=[Vol(1),Vol(2),Vol(3),Vol(4),Vol(1)+Vol(2)+Vol(3)+Vol(4)]';
                                    peso=[Pes(1),Pes(2),Pes(3),Pes(4),PesoRes]';
                                    brazo=[B(1),B(2),B(3),B(4)," "]';
                                    momento=[M(1),M(2),M(3),M(4),Mresistente]';
                                    Tabla= FormalTable(["SECCIÓN","VOLUMEN","PESO","BRAZO","MOMENTO"],[seccion,volumen,peso,brazo,momento]);
                                    Tabla.Style = {OuterMargin("0pt", "0pt","50pt","5pt")};
                                    Tabla.Header.Style{end+1} = BackgroundColor("silver");
                                    Tabla.Width = "500pt";
                                    Tabla.Border="solid";
                                    Tabla.ColSep="solid";
                                    Tabla.RowSep="solid";
                                    Tabla.HAlign="center";
                                    Tabla.TableEntriesHAlign="center";
                                    
                                end
                                if Cant==1
                                    image2 = mlreportgen.report.FormalImage();
                                    image2.Image = which('reportetriza.png');
                                    image2.Caption = 'Partes Resistente(Muro Canteliver)';
                                    image2.Height = '3in';
                                    image2.Width = '4in'; 

                                    seccion=["A","B","C","D"," "]';
                                    volumen=[Vol(1),Vol(2),Vol(3),Vol(4),Vol(1)+Vol(2)+Vol(3)+Vol(4)]';
                                    peso=[Pes(1),Pes(2),Pes(3),Pes(4),PesoRes]';
                                    brazo=[B(1),B(2),B(3),B(4)," "]';
                                    momento=[M(1),M(2),M(3),M(4),Mresistente]';
                                    Tabla= FormalTable(["SECCIÓN","VOLUMEN","PESO","BRAZO","MOMENTO"],[seccion,volumen,peso,brazo,momento]);
                                    Tabla.Style = {OuterMargin("0pt", "0pt","50pt","5pt")};
                                    Tabla.Header.Style{end+1} = BackgroundColor("silver");
                                    Tabla.Width = "500pt";
                                    Tabla.Border="solid";
                                    Tabla.ColSep="solid";
                                    Tabla.RowSep="solid";
                                    Tabla.HAlign="center";
                                    Tabla.TableEntriesHAlign="center";
                                
                                end


                                            
                                                             
                                add(s,image2);
                               
%                                 plot1 = Image("reportetriza.png");
%                                 plot1.Width = "4in";
%                                 plot1.Height = "3in";
%                                 plot1.Style = {OuterMargin("150pt", "0pt","0pt","5pt")};
%                                 append(s,plot1);
                                
                               

                                
                            
                               
                               
                                
                                
                                add(s,Tabla)
                                add(chap2,s)
                                add(rpt,chap2)
                                %%%%%%%Capitulo3
                                chap3 = Chapter('Factores Actuantes');
                                if k==2
                                    s3=Section("Condición Activa de Rankine");
                                     condicion=["Ka"]';
                                end
                                if k==3
                                    s3=Section("Condición Reposo de Rankine");
                                     condicion=["Ko"]';
                                end
                                


%                                 plot3 = Image("kacou.png");
%                                 plot3.Style = {OuterMargin("100pt", "0pt","0pt","5pt")};
%                                 append(s3,plot3);
                    
                                
                                valor=[ka]';
                               
                                
                                Tabla3= FormalTable(["CONDICIÓN", ""],[condicion,valor]);
                                Tabla3.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla3.Style = {OuterMargin("0pt", "0pt","30pt","30pt")};
                                Tabla3.Width = "200pt";
                                Tabla3.Border="solid";
                                Tabla3.ColSep="solid";
                                Tabla3.RowSep="solid";
                                Tabla3.HAlign="center";
                                Tabla3.TableEntriesHAlign="center";
                                Tabla3.Style= [Tabla3.Style {NumberFormat("%.2f")}];
                                add(s3,Tabla3);
                    
                                %%%%TABLA FUERZAS
                                s32=Section("Fuerzas");
                                secc=[1,2,3,"agua"]';
                                fue=[F,pesoagua]';
                                mom=[Maa,magua]';
                                fue=round(fue,3);
                                mom=round(mom,3);
                                
                                Tabla31= FormalTable(["SECCIÓN"," Fuerzas ","Momento "],[secc,fue,mom]);
                                Tabla31.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla31.Style = {OuterMargin("0pt", "0pt","40pt","40pt")};
                                Tabla31.Width = "400pt";
                                Tabla31.Border="solid";
                                Tabla31.ColSep="solid";
                                Tabla31.RowSep="solid";
                                Tabla31.HAlign="center";
                                Tabla31.TableEntriesHAlign="center";
                                Tabla31.Style= [Tabla31.Style {NumberFormat("%.2f")}];
                                add(s32,Tabla31);
%                                 %%%%FUERZAS DESCOMPUESTAS
%                                 secc1=["Pav (Pa*sen(δ))","Pah (Pa*cos(δ))"]';
%                                 fue=[pvertical,phor]';
%                                 bra=[Bvertical,Ba]';
%                                 mom=[Mvertical,Ma(1)]';
%                                 
%                                 Tabla32= FormalTable(["SECCIÓN"," Fuerzas ","Brazos","Momento "],[secc1,fue,bra,mom]);
%                                 Tabla32.Header.Style{end+1} = BackgroundColor("silver");
%                                 Tabla32.Style = {OuterMargin("0pt", "0pt","40pt","40pt")};
%                                 Tabla32.Width = "500pt";
%                                 Tabla32.Border="solid";
%                                 Tabla32.ColSep="solid";
%                                 Tabla32.RowSep="solid";
%                                 Tabla32.HAlign="center";
%                                 Tabla32.TableEntriesHAlign="center";
%                                 Tabla32.Style= [Tabla32.Style {NumberFormat("%.2f")}];
%                                 add(s32,Tabla32);
                               
                    
                                add(chap3,s3);
                                add(chap3,s32)
                                add(rpt,chap3);
                     %%%%%%%Capitulo5
                    
                    
                                chap5 = Chapter('Factores de seguridad estaticos');
                                s5=Section("Factor de seguridad al desplazamiento");
                                image5 = mlreportgen.report.FormalImage();
                                image5.Image = which('fsdes.png');
                                image5.Caption = 'Ecucaion para el Factor de seguridad al deslizamiento ';
                                add(s5,image5);
%                                 plot5 = Image("fsdes.png");
%                                 plot5.Style = {OuterMargin("100pt", "0pt","0pt","5pt")};
%                                 append(s5,plot5);
                                s52=Section("Factor de seguridad al volcamiento");
                                image52 = mlreportgen.report.FormalImage();
                                image52.Image = which('fsvol.png');
                                image5.Caption = 'Ecucaion para el Factor de seguridad al volcamiento ';
                                add(s52,image52);

%                                 plot52 = Image("fsvol.png");
%                                 plot52.Style = {OuterMargin("100pt", "0pt","0pt","5pt")};
%                                 append(s52,plot52);
%                                 
                                s53=Section("Resultados");
                                %TABLARESULTADOS
                                if So==0
                                    fsde=[FSD]';
                                    fsv1=[FSV]';
                                end
                                

                                if So>0
                                    %antes
                                    if aa==2
                                        fsde=[FSD1]';
                                        fsv1=[FSV1]';
                                    end
                                    %despues
                                    if aa==3
                                        fsde=[FSD2]';
                                        fsv1=[FSV2]';
                                    end

                                    %ambas
                                    if aa==4
                                        fsde=[FSD3]';
                                        fsv1=[FSV3]';
                                    end

                                end
                                
                                
                               
                                Tabla5= FormalTable(["F.S DES","F.S.1 VOL"],[fsde,fsv1]);
                                Tabla5.Style = {OuterMargin("0pt", "0pt","50pt","5pt")};
                                Tabla5.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla5.Width = "400pt";
                                Tabla5.Border="solid";
                                Tabla5.ColSep="solid";
                                Tabla5.RowSep="solid";
                                Tabla5.HAlign="center";
                                Tabla5.TableEntriesHAlign="center";
                                Tabla5.Style= [Tabla5.Style {NumberFormat("%.2f")}];
                                add(s53,Tabla5);
                    
                    
                                add(chap5,s5);
                                add(chap5,s52);
                                add(chap5,s53);
                                add(rpt,chap5);

                                % Close the report (required)
                                close(rpt)
                                % Display the report (optional)
                                rptview(rpt);
            
            
                             case 'NO'

                        end
            

            end

     
                                


end