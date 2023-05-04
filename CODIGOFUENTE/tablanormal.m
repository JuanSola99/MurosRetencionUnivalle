function [D,Mmax,Zmax,Ss]=tablanormal(espesor,phi,PesoE,n,NF,rr1,rr2,eper)
     if ismcc || isdeployed
        % Make sure DOM is compilable
        makeDOMCompilable();
     end
    import mlreportgen.dom.*
    import mlreportgen.report.*
  
    %%% METODO DE BLUM
    h=espesor;
    g=PesoE; 
    phii=phi;
    nf=NF;
    
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

     %%%%INICIO METODO DE BLUM 1 estrato uniforme
     re=0.6;
     Kpm=(tand(45+(phii)/2)^2)*re;
     HT=espesor;
     to=HT/((Kpm/ka)^(1/3)-1);
     dto=0.2*to;
     D=to+dto;%solo prof hincado
     Dtotal=D+HT;
     assignin('base','KPM',Kpm);
     assignin('base','to',to);
      %ESTA SECO
      if nf==0
             EF=g;
             assignin('base','EF',EF);
            
      end
      %ESTA SATURADO

      if nf==1
             EF=g-un;
             assignin('base','EF',EF);
            
      end

      

     %EMPUJE ACTIVO
     PA=ka*g*(HT+to)^2/2;

     %PASIVO
     PP=Kpm*g*(to)^2/2;

     %FUERZA C
     C=PP-PA;

     %MOMENTO
     Mmax=(1/6)*g*Kpm*((HT^3)/(sqrt(Kpm/ka)-1)^2);
     %Prof mom max si se le suma h es para que quede desde la superficie
     Zmax=(HT/(sqrt(Kpm/ka)-1));

     %Prof mom max si se le suma h es para que quede desde la superficie
     Zmax2=(HT/(sqrt(Kpm/ka)-1)+HT);



    
     %MMAX dependiente de zmax
     Mmaxz=(1/6)*g*(HT+Zmax)^3*ka-((1/6)*g*Zmax^3*kp);

     

    
     

     %Seccion tanbla

     Ss=Mmax/eper;
     assignin('base','Ss',Ss)


    ;
     assignin('base','Zmax',Zmax2)

     assignin('base','C',C);
     assignin('base','Mmax',Mmax);
     assignin('base','Mmaxz',Mmaxz);

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
                                titlepg.Title = 'ANÁLISIS GEOTÉCNICO DE TABLESTACAS MÉTODO DE BLUM';
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
                                seccion1=["TIPO DE TABLESTACA","ALTURA TABLESTACA ANTES DEL HINCADO"," COEFICIENTE DE SEGURIDAD (ϒe)","ESFUERZO PERMISIBLE DEL MATERIAL"]';
                                volumen1=["TABLESTACA EN VOLADIZO", h ,"0,6", eper]';
                                
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
                                if nf==1
                                    seccions1=["PESO ESPECIFICO (suelo saturado)","ANGULO DE FRICCION (ϕ)","COHESIÓN"]';
                                end
                                if nf==0
                                    seccions1=["PESO ESPECIFICO (suelo seco)","ANGULO DE FRICCION (ϕ)","COHESIÓN"]';
                                end
                                
                                volumens1=[PesoE,phi,0]';
                               
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
                                image3.Image = which('esfuer.png');
                                image3.Caption = 'Pantalla en voladizo. [Blum 1951] ';

                                 add(s3,image3);

                                image31 = mlreportgen.report.FormalImage();
                                image31.Image = which('kpm.png');
                                image31.Caption = 'Coeficientes de condición del terreno';
                                

                                 add(s3,image31);
                    
                                condicion=["Ka","Kp","Kpm"]';
                                valor=[ka,kp,Kpm]';
                               
                                
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

                                s32=Section("Profundiad de Hincado");

                                image32 = mlreportgen.report.FormalImage();
                                image32.Image = which('to.png');
                                image32.Caption = 'Punto de giro de la tablestaca ';

                                 add(s32,image32);

                                image33 = mlreportgen.report.FormalImage();
                                image33.Image = which('dto.png');
                                image33.Caption = 'Punto de giro de la tablestaca ';

                                 add(s32,image33);

                                condicion=["To","ΔTo","D"]';
                                valor=[to,dto,D]';
                               
                                
                                Tabla32= FormalTable(["PROFUNDIDADES","VALOR"],[condicion,valor]);
                                Tabla32.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla32.Style = {OuterMargin("0pt", "0pt","30pt","30pt")};
                                Tabla32.Width = "200pt";
                                Tabla32.Border="solid";
                                Tabla32.ColSep="solid";
                                Tabla32.RowSep="solid";
                                Tabla32.HAlign="center";
                                Tabla32.TableEntriesHAlign="center";
                                add(s32,Tabla32);

                                add(chap3,s3)
                                add(chap3,s32)


                                add(rpt,chap3);


                                %%%% Chap 4

                                 %%%%TABLA FUERZAS
                                chap4= Chapter('Factores Actuantes');
                                s4=Section("Emupujes");
                                image4 = mlreportgen.report.FormalImage();
                                image4.Image = which('empuje.png');
                                image4.Caption = 'Formula para el empuje activo y pasivo';
                                

                                 add(s4,image4);


                                secc=["Empuje activo","Empuje pasivo ","Fuerza C"]';
                                fue=[PA,PP,C]';
                                
                                
                                Tabla31= FormalTable(["EMPUJE","VALOR"],[secc,fue]);
                                Tabla31.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla31.Style = {OuterMargin("0pt", "0pt","40pt","40pt")};
                                Tabla31.Width = "400pt";
                                Tabla31.Border="solid";
                                Tabla31.ColSep="solid";
                                Tabla31.RowSep="solid";
                                Tabla31.HAlign="center";
                                Tabla31.TableEntriesHAlign="center";
                                add(s4,Tabla31);
                                add(chap4,s4)

                                add(rpt,chap4);

                                %%%%%%%%%%%%%
                                chap5= Chapter('Momento maximo y donde ocurre');
                                s5=Section("Momento Maximo");
                                image5 = mlreportgen.report.FormalImage();
                                image5.Image = which('mmax.png');
                                image5.Caption = 'Formula para el Momento maximo y el punto donde ocurre desde la superficie de terrenoclc';
                                

                                 add(s5,image5);


                                secc=["Mmax","Zmax"]';
                                fue=[Mmax,Zmax2]';
                                
                                
                                Tabla5= FormalTable([" ","VALOR"],[secc,fue]);
                                Tabla5.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla5.Style = {OuterMargin("0pt", "0pt","40pt","40pt")};
                                Tabla5.Width = "400pt";
                                Tabla5.Border="solid";
                                Tabla5.ColSep="solid";
                                Tabla5.RowSep="solid";
                                Tabla5.HAlign="center";
                                Tabla5.TableEntriesHAlign="center";
                                add(s5,Tabla5);

                                add(chap5,s5)

                                add(rpt,chap5);


                                 % Close the report (required)
                                close(rpt)
                                % Display the report (optional)
                                rptview(rpt);





                            case 'NO'
                                
                                
                            
                        end









end