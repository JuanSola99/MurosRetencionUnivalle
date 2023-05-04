function [Dd,FT,MMAX,Ff,Ss,PS,LP]=tablaanclado(espesor,phi,PesoE,n,NF,AA,rr1,rr2,l1,l2,FS,eper)
     if ismcc || isdeployed
        % Make sure DOM is compilable
        makeDOMCompilable();
     end
    import mlreportgen.dom.*
    import mlreportgen.report.*
  
    %%% METODO DE braja das Método de apoyo empotrado en tierra para penetración en suelo arenoso
    assignin('base','l1',l1);
    assignin('base','l2',l2);
    h=espesor;
    sumh=sum(h);
    g=PesoE; 
    phii=phi;
    nf=NF;
    assignin('base','rrr1',rr1);
    %%%%%%%
    ka=zeros(n,1);
    if rr1==1
        un=9.81;
        
    end
    if rr2==1
        un=62.4;
    end

    assignin('base','un',un);
        
    for j=1:n
        ka(j)=(1-sind(phii(j)))/(1+sind(phii(j)));
       

    end
     assignin('base','ka',ka);

    kp=zeros(n,1);
    for j=1:n
        kp(j)=(1+sind(phii(j)))/(1-sind(phii(j)));
        

    end
     assignin('base','kp',kp);
    
    
    
    %%%%%%%VERTICAL
    
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
            sv(kk2)=h(kk2)*(g(kk2)-un);
        end
    end
    if nf==0
        %ESTA SECO
        for kk=1:n
            sv(kk)=h(kk)*g(kk);
        end
        
    end
    
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
        end
    end
    
     assignin('base','Sa',sa);
     assignin('base','SB',sb);
     assignin('base','SS',SS);

     

     %%%%%%%%%no hay nivel freatico
     if n==1
         if nf==0
             EF=g;
                        

         end

         if nf==1
             warndlg('Para este caso donde el # de estratos es 1 no se tiene encuenta un nivel freatico, asumir nf=0','Observación')
             return

         end
          kpa=kp-ka; 
         
         assignin('base','EF',EF);
          assignin('base','kpa',kpa);

         if phii==30
             L5=0.08*(h);

         end
         if phii==35

             L5=0.03*(h1);

         end
         if phii==40
             L5=0*h;

         end
         assignin('base','L5',L5);
          %%%PRESION ACTIVA NETA DEBAJO DE LA LINEA DE DRAGADO
             PEF=SS(2)-(EF*kpa*L5);
             assignin('base','PEF',PEF);

              %%%L'
             LP=L5+(h-l1);
             assignin('base','LP',LP);

             %%%%esfuerzos por ancalaje 
             l12=l1+l2;
             
             s1=l1*ka(1)*g(1);
             s2=l12*ka(1)*g(1);
             
             S2=(EF*h)*ka;
             S2=S2(1);

             assignin('base','s1',s1);
             assignin('base','s2',s2);
             assignin('base','S2',S2);

             %%%PRESION DIAGRAMA DE PRESION 1 ESTRATOS
             h2=h-l12;

             w1=(0.5*(s1+s2)*l2)+0.5*h2*(s2+S2)+(0.5*L5*(S2+PEF));
             assignin('base','w1',w1);




             %%%MOMENTO
             
             MMAX=w1*LP/8;
             MMAX=MMAX(1);
             assignin('base','MMAX',MMAX);


             %PPRIMA
             %brazo respecto a o'
             b1=((2/3*l12)-l1);
             b2=(l2+h2/2);
             b3=(l2+2/3*h2);
             b4=L5*(l2+h2+L5/2);
             assignin('base','b1',b1);
             %momento del área ACDJI
             p1=0.5*s2*l12;
             p2=s2*h2;
             p3=0.5*h2*(S2-s2);
             p4=0.5*(S2+PEF);
             assignin('base','p1',p1);

             %
             M1=p1*b1;
             M2=p2*b2;
             M3=p3*b3;
             M4=p4*b4;

             assignin('base','M1',M1);
             MACD=M1+M2+M3+M4;
             assignin('base','MACD',MACD);
             PP=(1/LP)*MACD;
             PP=PP(1);
             assignin('base','PP',PP);

             %PROF HINCADO
             divi=6*PP/(EF*kpa);
             assignin('base','Divi',divi);
             Dd=L5+(1.2*divi^0.5);
             assignin('base','Dd',Dd);

             %TOMANDO MOMENTO RESPECTO A I
             bb1=L5+h2+(l12/2);
             bb2=L5+(h2/2);%parte cuadrada
             bb3=L5+(h2)/3;%parte trinagular
             bb4=L5/2;
             
             pp1=0.5*s2*l12;
             pp2=s2*h2;
             pp3=0.5*h2*(S2-s2);
             pp4=0.5*(PEF+S2)*L5;

             assignin('base','pp1',pp1);
             assignin('base','pp3',pp3);
             assignin('base','pp4',pp4);


             MM1=pp1*bb1;
             MM2=pp2*bb2;
             MM3=pp3*bb3;
             MM4=pp4*bb4;

             MI=MM1+MM2+MM3+MM4;

             FT=(1/LP)*MI;
             FT=FT(1);
             assignin('base','FT',FT);

             Ff=FT/cosd(AA);


         
     end

     %%%%INICIO METODO DE Braja Das 1 estrato con nf=2 

     if n==2
         

         if nf==2
             EF=g(2)-un;
             assignin('base','EF',EF);
             kpa=(kp(2)-ka(2));
             assignin('base','kpa',kpa);
         end
        

             if phii==30
                 L5=0.08*(h(1)+h(2));

             end
             if phii==35
                 L5=0.03*(h(1)+h(2));

             end
             if phii==40
                 L5=0*(h(1)+h(2));

             end
             assignin('base','L5',L5);

             %%%PRESION ACTIVA NETA DEBAJO DE LA LINEA DE DRAGADO
             PEF=SS(3)-(EF*kpa*L5);
             assignin('base','PEF',PEF);

             %%%L'
             LP=l2+h(2)+L5;
             assignin('base','LP',LP);

             %%%%esfuerzos por ancalaje 
             l12=l1+l2;
             
             s1=l1*ka(1)*g(1);
             s2=l12*ka(1)*g(1);
             
             S2=(g(1)*l12+EF*h(2))*ka;
             S2=S2(1);

             assignin('base','s1',s1);
             assignin('base','s2',s2);
             assignin('base','S2',S2);

             %%%PRESION DIAGRAMA DE PRESION 2 ESTRATOS
             h2=h(2);

             w1=(0.5*(s1+s2)*l2)+0.5*h2*(s2+S2)+(0.5*L5*(S2+PEF));
             assignin('base','w1',w1);




             %%%MOMENTO
             
             MMAX=w1*LP/8;
             MMAX=MMAX(1);
             assignin('base','MMAX',MMAX);


             %PPRIMA
             %brazo respecto a o'
             b1=((2/3*l12)-l1);
             b2=(l2+h(2)/2);
             b3=(l2+2/3*h(2));
             b4=L5*(l2+h(2)+L5/2);
             assignin('base','b1',b1);
             %momento del área ACDJI
             p1=0.5*s2*l12;
             p2=s2*h(2);
             p3=0.5*h(2)*(S2-s2);
             p4=0.5*(S2+PEF);
             assignin('base','p1',p1);

             %
             M1=p1*b1;
             M2=p2*b2;
             M3=p3*b3;
             M4=p4*b4;

             assignin('base','M1',M1);
             MACD=M1+M2+M3+M4;
             assignin('base','MACD',MACD);
             PP=(1/LP)*MACD;
             PP=PP(1);
             assignin('base','PP',PP);

             %PROF HINCADO
             divi=6*PP/(EF*kpa);
             assignin('base','Divi',divi);
             Dd=L5+(1.2*divi^0.5);
             assignin('base','Dd',Dd);

             %TOMANDO MOMENTO RESPECTO A I
             bb1=L5+h(2)+(l12/2);
             bb2=L5+(h(2)/2);%parte cuadrada
             bb3=L5+(h(2))/3;%parte trinagular
             bb4=L5/2;
             
             pp1=0.5*s2*l12;
             pp2=s2*h(2);
             pp3=0.5*h(2)*(S2-s2);
             pp4=0.5*(PEF+S2)*L5;

             assignin('base','pp1',pp1);
             assignin('base','pp3',pp3);
             assignin('base','pp4',pp4);


             MM1=pp1*bb1;
             MM2=pp2*bb2;
             MM3=pp3*bb3;
             MM4=pp4*bb4;

             MI=MM1+MM2+MM3+MM4;

             FT=(1/LP)*MI;
             FT=FT(1);
             assignin('base','FT',FT);

             Ff=FT/cosd(AA);

            


          



     end

      PS=Dd*FS;

     Ss=MMAX/eper;
     assignin('base','Ss',Ss)

     %%%%%%%%%%%%%%%GENERAR REPORTE

         reporte=questdlg('Desea generar reporte?', ...
	                        'Generar reporte', ...
	                        'SI',' NO  ','  CANCELAR  ','');
                        % Handle response
                        switch reporte
                            case 'SI'   
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
                                titlepg.Title = 'ANÁLISIS GEOTÉCNICO DE TABLESTACAS MÉTODO DE BRAJA DAS. (basado en Cornfield, 1975)';
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
                               
                                %TABLA MURO
                                s1=Section("CARACTERISTICAS MURO");
                                seccion1=["TIPO DE TABLESTACA","ALTURA TABLESTACA ANTES DEL HINCADO","ESFUERZO PERMISIBLE DEL MATERIAL"]';
                                volumen1=["TABLESTACA CON ANCLAJE",sumh,eper]';
                                
                                Tabla1= FormalTable(["PROPIEDAD","VALOR"],[seccion1,volumen1]);
                                Tabla1.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla1.Width = "500pt";
                                Tabla1.Border="solid";
                                Tabla1.ColSep="solid";
                                Tabla1.RowSep="solid";
                                Tabla1.HAlign="left";
                                Tabla1.TableEntriesHAlign="center";

                                add(s1,Tabla1);
                                add(chap,s1);

                                
                                %TABLA SUELO RELLENO
                                ss1=Section("CARACTERISTICAS SUELO DE RELLENO");
                                if n==1
                                   
                                    
                                    if nf==0
                                    seccions1=["PESO ESPECIFICO (suelo seco)","ANGULO DE FRICCION (ϕ)","COHESIÓN"]';
                                    volumens1=[PesoE,phi,0]';
                                    end
                                end

                                 if n==2
                                    if nf==2
                                    seccions1=["PESO ESPECIFICO (estrato 1 )","PESO ESPECIFICO (estrato 2 saturado)","ANGULO DE FRICCION (ϕ1)","ANGULO DE FRICCION (ϕ2)","COHESIÓN"]';
                                    volumens1=[PesoE(1),PesoE(2),phi(1),phi(2),0]';
                                    end
                                    %para n==2 solo se admite nf==2
                                end
                                
                                
                                
                               
                                 %Base1,Base2,Ha,D,B,alpha,phib,concreto,espesor,phi,PesoE,Cohesion,n,k,espe,pe,ph,R,TS,R1,R2,kh,kv,So
                                
                                Tablas1= FormalTable(["PROPIEDAD","VALOR"],[seccions1,volumens1]);
                                Tablas1.Header.Style{end+1} = BackgroundColor("silver");
                                Tablas1.Width = "500pt";
                                Tablas1.Border="solid";
                                Tablas1.ColSep="solid";
                                Tablas1.RowSep="solid";
                                Tablas1.HAlign="left";
                                Tablas1.TableEntriesHAlign="center";

                                add(ss1,Tablas1);
                                add(chap,ss1);                                                         
                                                                                
                               
                                                                                  
                                
                                add(rpt,chap);

                                %%%%%%%Capitulo3
                                chap3 = Chapter('Factores Actuantes');
                                s3=Section("Condición Activa y Pasiva sobre la tablestaca");

                                image3 = mlreportgen.report.FormalImage();
                                image3.Image = which('fbraja.png');
                                image3.Caption = 'Diagrama de esfuerzos ';

                                 add(s3,image3);

                                image31 = mlreportgen.report.FormalImage();
                                image31.Image = which('kpb.png');
                                image31.Caption = 'Coeficientes de condición del terreno';
                                

                                 add(s3,image31);
                    
                                condicion=["Ka","Kp"]';
                                valor=[ka(1),kp(1)]';
                               
                                
                                Tabla3= FormalTable(["CONDICIÓN", " "],[condicion,valor]);
                                Tabla3.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla3.Style = {OuterMargin("0pt", "0pt","30pt","30pt")};
                                Tabla3.Width = "200pt";
                                Tabla3.Border="solid";
                                Tabla3.ColSep="solid";
                                Tabla3.RowSep="solid";
                                Tabla3.HAlign="center";
                                Tabla3.TableEntriesHAlign="center";
                                add(s3,Tabla3);

                                  %%%%secccion32

                                s32=Section("Presion efectiva neta");
                                image322 = mlreportgen.report.FormalImage();
                                image322.Image = which('l5.png');
                                image322.Caption = 'Valor de L5';
                                

                                 add(s32,image322);

                                image32 = mlreportgen.report.FormalImage();
                                image32.Image = which('pef.png');
                                image32.Caption = 'Presión efectiva neta a una profundidad L5 ';

                              
                                 add(s32,image32);

                                  %%%%secccion33

                                s33=Section("Momento Maximo");
                                image3 = mlreportgen.report.FormalImage();
                                image3.Image = which('l5.png');
                                image3.Caption = 'Valor de L5';
                                

                                 add(s33,image3);

                                image33 = mlreportgen.report.FormalImage();
                                image33.Image = which('pef.png');
                                image33.Caption = 'Presión efectiva neta a una profundidad L5 ';

                               add(s33,image33);

                                image333 = mlreportgen.report.FormalImage();
                                image333.Image = which('pesob.png');
                                image333.Caption = 'Peso';
                                add(s33,image333);

                                image3333 = mlreportgen.report.FormalImage();
                                image3333.Image = which('mmaxb.png');
                                image3333.Caption = 'Momento maximo';
                                add(s33,image3333);

                                secc=["Empuje Efectivo(PEF)","Presión(W)","Momento Maximo(Mmax)"]';
                                fue=[PEF,w1,MMAX]';
                                
                                
                                Tabla31= FormalTable(["EMPUJE","VALOR"],[secc,fue]);
                                Tabla31.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla31.Style = {OuterMargin("0pt", "0pt","40pt","40pt")};
                                Tabla31.Width = "400pt";
                                Tabla31.Border="solid";
                                Tabla31.ColSep="solid";
                                Tabla31.RowSep="solid";
                                Tabla31.HAlign="center";
                                Tabla31.TableEntriesHAlign="center";
                                add(s33,Tabla31);



                                add(chap3,s3)
                                add(chap3,s32)
                                add(chap3,s33)    

                                add(rpt,chap3);

                                %%%% Chap 4

                                chap4= Chapter('PRESIÓN NETA, PROFUNDIAD DE HINCADO Y FUERZA DE ANCLAJE');
                                s4=Section("Presión neta");
                                image4 = mlreportgen.report.FormalImage();
                                image4.Image = which('pneta.png');
                                image4.Caption = 'Diagrama de presón neta';
                                

                                 add(s4,image4);

                                image41 = mlreportgen.report.FormalImage();
                                image41.Image = which('pp.png');
                                image41.Caption = 'Presión Neta';
                                

                                 add(s4,image41);

                                s41=Section("Profundidad de anclaje");
                                image42 = mlreportgen.report.FormalImage();
                                image42.Image = which('dbraja.png');
                                image42.Caption = 'Profundiad de anclaje D';

                                add(s41,image42);

                                %%%%

                                s42=Section("FUERZA DE ANCLAJE");
                                image43 = mlreportgen.report.FormalImage();
                                image43.Image = which('fbra.png');
                                image43.Caption = 'F';

                                add(s42,image43);

                                secc=["Presión Neta(P')","Profundiad de Hincado(D)","Fuerza de Anclaje(F)"]';
                                fue=[PP,Dd,FT]';
                                
                                
                                Tabla4= FormalTable([" ","VALOR"],[secc,fue]);
                                Tabla4.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla4.Style = {OuterMargin("0pt", "0pt","40pt","40pt")};
                                Tabla4.Width = "400pt";
                                Tabla4.Border="solid";
                                Tabla4.ColSep="solid";
                                Tabla4.RowSep="solid";
                                Tabla4.HAlign="center";
                                Tabla4.TableEntriesHAlign="center";
                                add(s42,Tabla4);
                                add(chap4,s4)

                                add(chap4,s41);
                                add(chap4,s42);

                                add(rpt,chap4);


                                 % Close the report (required)
                                close(rpt)
                                % Display the report (optional)
                                rptview(rpt);

                               


                            case 'NO'
                                
                                
                            
                        end






end




    







