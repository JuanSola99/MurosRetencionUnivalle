function [FSD,FSVR,FSVA]= Mtriangular(Base1,Base2,Ha,D,B,alpha,phib,concreto,espesor,phi,PesoE,Cohesion,n,k,espe,pe,ph,R,TS,R1,R2,Carga)

     if ismcc || isdeployed
        % Make sure DOM is compilable
        makeDOMCompilable();
     end
    import mlreportgen.dom.*
    import mlreportgen.report.*

    h=espesor;
    g=PesoE;
    c=Cohesion;   
    phii=phi;
    a=0;
    %Si hay suelo posterior
    assignin('base','R',R);
    assignin('base','alpha',alpha);
    if R==1
        a=0;
        nume=sind(B-ph)^2;
        p11=(sind(B)^2)*sind(B+D);
        p22=(1-sqrt(sind(ph+D)*sind(ph+a)/(sind(B+D)*sind(a+B))))^2;
        deno=p11*p22;
        kp=nume/deno;
        assignin('base','kp',kp);
        %Esfuerzo vertical
        ssv=espe*pe;
        assignin('base','espe',espe);
        assignin('base','pe',pe);
        %Esfuerzo horizontal
        ssh=ssv*kp;
        assignin('base','ssh',ssh);
        %Fuerza
        Fp=ssh*espe/2;
        Fps=Fp*sind(D);
        Fpc=Fp*cosd(D);
        %Brazo
        brap=espe/3;
        bras=Base2*2/3;
        brac=Ha/3;
        assignin('base','Fp',Fp);
        assignin('base','Fpc',Fpc);
        assignin('base','Fps',Fps);
        assignin('base','brap',brap);
        assignin('base','bras',bras);
        assignin('base','brac',brac);
        %Momento
        Mp=Fp*brap;
        assignin('base','Mp',Mp);
        
        Mps=Fps*bras;
        Mpc=Fpc*brac;

    end
    
    %Ka coulomb 
    %hallar ka o ko 
    if k==2    %cuando es activo
        num=sind(B+phii)^2;
        p1=(sind(B)^2)*sind(B-D);
        p2=(1+sqrt(sind(phii+D)*sind(phii-alpha)/(sind(B-D)*sind(alpha+B))))^2;
        den=p1*p2;
        ka=num/den;
        
        
    end
    if k==3 %cuando es pasivo
        nume=sind(B-ph)^2;
        p11=(sind(B)^2)*sind(B+D);
        p22=(1-sqrt(sind(ph+D)*sind(ph+a)/(sind(B+D)*sind(a+B))))^2;
        deno=p11*p22;
        ka=nume/deno;
        
    end
     assignin('base','Ka',ka);
   
    
    %ACTUANTE
    if n==1
        if TS==3
            if R1==1
                Pactuante=1/2*(g-9.81)*h^2*ka;
            end
            if R2==1
                Pactuante=1/2*(g-62.4)*h^2*ka;
            end
             %%%%%%%%%%%%%%%%%%%CARGA
            if Carga==0
                Pcarga=0;
                Pactuante=Pactuante+Pcarga;
                
            end
            if Carga>0
                factor=sind(B)^2/abs(sind(B+alpha));
                Pcarga=ka*Carga*factor*h;
                Pactuante=Pactuante+Pcarga;
                assignin('base','factor',factor);
               
            end
             assignin('base','Pactuante',Pactuante);
             assignin('base','Pcarga',Pcarga);
            %%%%%%%%%%%%%%%5%%
            pvertical=Pactuante*sind(D);
            phor=Pactuante*cosd(D);


           
            if R1==1
                Pesoagua=1/2*9.81*h^2;
            end
            if R2==1
                Pesoagua=1/2*62.4*h^2;
            end
            
            %Brazo
            Ba=h/3;
            Bagua=h/3;
            Bvertical=Base1+Base2;
            %Momento
            Ma(1)=phor*Ba;
            Ma(2)=Pesoagua*Bagua;
            Mact=sum(Ma);
            Mvertical=Bvertical*pvertical;
            
            
            
        end
        if TS==2
            Pactuante=1/2*(g)*h^2*ka;
            %%%%%%%%%%%%%% Carga

            if Carga==0
                Pcarga=0;
                Pactuante=Pactuante+Pcarga;
                
            end
            if Carga>0
                factor=sind(B)^2/abs(sind(B+alpha));
                Pcarga=ka*Carga*factor*h;
                Pactuante=Pactuante+Pcarga;
                assignin('base','factor',factor);
               
            end
            %%%%%%%%%%%%%%%%%%
            pvertical=Pactuante*sind(D);
            phor=Pactuante*cosd(D);
            Pesoagua=0;
            
             assignin('base','Pactuante',Pactuante);
             assignin('base','Pcarga',Pcarga);
            %Brazo
            Ba=h/3;
            Bagua=h/3;
            Bvertical=Base1+Base2;
            %Momento
            Ma(1)=phor*Ba;
            Ma(2)=Pesoagua*Bagua;
            Mact=sum(Ma);
            Mvertical=Bvertical*pvertical;
        end
        assignin('base','Mact',Mact);
        assignin('base','Bvertical',Bvertical);
         assignin('base','pvertical',pvertical);
         assignin('base','phor',phor);
         assignin('base','Mvertical',Mvertical);
                
    end
    
    %RESISTENTE
    %Volumen
    Vol(1)=Base1*Ha*1;
    Vol(2)=Base2*Ha*1/2;
    %Peso
    Pes(1)=Vol(1)*concreto;%Kn
    Pes(2)=Vol(2)*concreto;%Kn
    %Brazos
    Br(1)=Base2+Base1/2;
    Br(2)=Base2*2/3; 
    %Momentos
    M(1)=Br(1)*Pes(1);
    M(2)=Br(2)*Pes(2);
    Mr=sum(M);
    
    
    PesoRes=sum(Pes);
    %FACTORES DE SEGURIDAD
    FSD=(((PesoRes+pvertical)*tand(phib))+(c*(Base1+Base2)))/(phor+ Pesoagua);
    FSVR=(Mr+Mvertical)/(Mact);
    FSVA=(Mr)/(Mact-Mvertical);
    
    if R==1
        FSD=((PesoRes+pvertical)*tand(phib))/(phor+ Pesoagua-Fpc);
        FSVR=(Mr+Mvertical)/(Mact);
        FSVA=Mr/(Mact-Mvertical);
    end
    assignin('base','FSD',FSD);
    assignin('base','FSVR',FSVR);
    assignin('base','FSVA',FSVA);
    

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
                               
                                %TABLA MURO
                                s1=Section("CARACTERISTICAS MURO");
                                seccion1=["ALTURA MURO","BASE DEL MURO (B1)","BASE DEL MURO (B2)","PESO ESPECIFICO (concreto)","ANGULO DE INCLINACIÓN DEL MURO (β)"]';
                                volumen1=[Ha,Base1,Base2,concreto,B]';
                                %Base1,Base2,Ha,D,B,alpha,phib,concreto,espesor,phi,PesoE,Cohesion,n,k,espe,pe,ph,R,TS,R1,R2,kh,kv,So
                                
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
                                if TS==3
                                    seccions1=["PESO ESPECIFICO (suelo saturado)","ANGULO DE FRICCION (ϕ1)","ANGULO DE FRICCION ENTRE EL SUELO Y EL MURO (δ)","INCLINACION DEL RELLENO (α)","COHESIÓN","SOBRE-CARGA"]';
                                end
                                if TS==2
                                    seccions1=["PESO ESPECIFICO (suelo seco)","ANGULO DE FRICCION (ϕ1)","ANGULO DE FRICCION ENTRE EL SUELO Y EL MURO (δ)","INCLINACION DEL RELLENO (α)","COHESIÓN","SOBRE-CARGA"]';
                                end
                                
                                volumens1=[PesoE,phi,D,alpha,0,Carga]';
                               
                                 %Base1,Base2,Ha,D,B,alpha,phib,concreto,espesor,phi,PesoE,Cohesion,n,k,espe,pe,ph,R,TS,R1,R2,kh,kv,So
                                
                                Tablas1= FormalTable(["PROPIEDAD","VALOR"],[seccions1,volumens1]);
                                Tablas1.Header.Style{end+1} = BackgroundColor("silver");
                                Tablas1.Width = "500pt";
                                Tablas1.Border="solid";
                                Tablas1.ColSep="solid";
                                Tablas1.RowSep="solid";
                                Tablas1.HAlign="left";
                                Tablas1.TableEntriesHAlign="center";

                                 

                                  %%%%%%%%%%%%%%%%%%%si hay suelo posterior
                                 if R==1

                                      seccionsp=["PESO ESPECIFICO (suelo posterior)","ANGULO DE FRICCION (ϕposterior)","ESPESOR"]';
                                      volumensp=[PesoE,ph,espe]';
                                      Tablasp= FormalTable(["PROPIEDAD SUELO POSTERIOR","VALOR"],[seccionsp,volumensp]);
                                      Tablasp.Header.Style{end+1} = BackgroundColor("silver");
                                      Tablasp.Width = "500pt";
                                      Tablasp.Border="solid";
                                      Tablasp.ColSep="solid";
                                      Tablasp.RowSep="solid";
                                      Tablasp.HAlign="left";
                                      Tablasp.TableEntriesHAlign="center";
                                      add(ss1,Tablasp);
                                      



                                 end

                    
                                %TABLA SUELO DE CIMENTACION
                                ssr1=Section("CARACTERISTICAS SUELO DE CIMENTACIÓN");
                                seccionsr1=["ANGULO DE FRICCION Base (ϕ2)","PESO ESPECIFICO (cimentación)","COHESIÓN"]';
                                volumensr1=[phib,PesoE,Cohesion]';
                               %phib,concreto,espesor,phi,PesoE,Cohesion,n,k,espe,pe,ph,R,TS,R1,R2,kh,kv,So
                                
                                
                                Tablasr1= FormalTable(["PROPIEDAD","VALOR"],[seccionsr1,volumensr1]);
                                Tablasr1.Header.Style{end+1} = BackgroundColor("silver");
                                Tablasr1.Width = "500pt";
                                Tablasr1.Border="solid";
                                Tablasr1.ColSep="solid";
                                Tablasr1.RowSep="solid";
                                Tablasr1.HAlign="left";
                                Tablasr1.TableEntriesHAlign="center";

                                add(ss1,Tablas1);
                                add(chap,ss1);

                                add(ssr1,Tablasr1);
                                
                                
                                add(chap,ssr1);                           
                               
                              
                                                          
                                
                                add(rpt,chap);
                                
                                

                               
                    
                                 %%%%%%%Capitulo2
                                chap2 = Chapter('Factores Resistentes');
                                s=Section("Partes resistentes");

                                image2 = mlreportgen.report.FormalImage();
                                image2.Image = which('reportetria.png');
                                image2.Caption = 'Partes Resistente(Muro Triangular)';
                                image2.Height = '3in';
                                image2.Width = '4in';                             
                                                             
                                add(s,image2);
%                                 plot1 = Image("reportetria.png");
%                                 plot1.Width = "4in";
%                                 plot1.Height = "3in";
%                                 plot1.Style = {OuterMargin("150pt", "0pt","0pt","5pt")};
%                                 append(s,plot1);

                                
                                
                                seccion=["A","B"," "]';
                                volumen=[Vol(1),Vol(2),Vol(2)+Vol(1)]';
                                peso=[Pes(1),Pes(2),Pes(2)+Pes(1)]';
                                brazo=[round(Br(1),3),round(Br(2),3)," "]';
                                momento=[M(1),M(2),Mr]';
                                
                                peso=round(peso,3);
                                
                                momento=round(momento,3);
                                Tabla= FormalTable(["SECCIÓN","VOLUMEN","PESO","BRAZO","MOMENTO"],[seccion,volumen,peso,brazo,momento]);
                                Tabla.Style = {OuterMargin("0pt", "0pt","50pt","5pt")};
                                Tabla.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla.Width = "500pt";
                                Tabla.Border="solid";
                                Tabla.ColSep="solid";
                                Tabla.RowSep="solid";
                                Tabla.HAlign="center";
                                Tabla.TableEntriesHAlign="center";

                                 %%%%%%%%%%%%%%%%%55posteriro

                                if R==1
                                     seccionp=["Posterior","Posterior horizontal"]';
                                     esfuerzop=[ssv,"",]';
                                     pesop=[Fp,Fpc]';
                                     brazop=[brap,""]';
                                     momentop=[Mp,"No se considera para el F.S"]';
                                     Tablap= FormalTable(["SECCIÓN","ESFUERZO","FUERZA","BRAZO","MOMENTO"],[seccionp,esfuerzop,pesop,brazop,momentop]);
                                     Tablap.Style = {OuterMargin("0pt", "0pt","50pt","5pt")};
                                     Tablap.Header.Style{end+1} = BackgroundColor("silver");
                                     Tablap.Width = "500pt";
                                     Tablap.Border="solid";
                                     Tablap.ColSep="solid";
                                     Tablap.RowSep="solid";
                                     Tablap.HAlign="center";
                                     Tablap.TableEntriesHAlign="center";
                                     add(s,Tablap);


                                end
                                
                                
                                
                                add(s,Tabla)
                                add(chap2,s)
                                add(rpt,chap2)
                                
                                
                                
                                
                                %%%%%%%Capitulo3
                                chap3 = Chapter('Factores Actuantes');
                                s3=Section("Condición Activa Coulomb");

                                image3 = mlreportgen.report.FormalImage();
                                image3.Image = which('kacou.png');
                                image3.Caption = 'Coeficiente Ka Coulomb ';
                                
%                                 plot3 = Image("kacou.png");
%                                 plot3.Style = {OuterMargin("100pt", "0pt","0pt","5pt")};
                                 add(s3,image3);
                    
                                condicion=["Ka"]';
                                valor=[ka]';
                                valor=round(valor,3);
                                
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
                    
                                %%%%TABLA FUERZAS
                                s32=Section("Fuerzas descompuestas");
                                secc=["Pa","Agua"]';
                                fue=[Pactuante,Pesoagua]';
                                mom=[" ",Ma(2)]';

                                fue=round(fue,3);
                                
                               
                                
                                Tabla31= FormalTable([" "," Fuerzas ","Momento "],[secc,fue,mom]);
                                Tabla31.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla31.Style = {OuterMargin("0pt", "0pt","40pt","40pt")};
                                Tabla31.Width = "400pt";
                                Tabla31.Border="solid";
                                Tabla31.ColSep="solid";
                                Tabla31.RowSep="solid";
                                Tabla31.HAlign="center";
                                Tabla31.TableEntriesHAlign="center";
                                add(s3,Tabla31);
                                %%%%FUERZAS DESCOMPUESTAS
                                secc1=["Pav (Pa*sen(δ))","Pah (Pa*cos(δ))"]';
                                fue=[pvertical,phor]';
                                bra=[Bvertical,Ba]';
                                mom=[Mvertical,Ma(1)]';
                                
                                fue=round(fue,3);
                                bra=round(bra,3);
                                mom=round(mom,3);

                                
                                Tabla32= FormalTable(["Componentes"," Fuerzas ","Brazos","Momento "],[secc1,fue,bra,mom]);
                                Tabla32.Header.Style{end+1} = BackgroundColor("silver");
                                Tabla32.Style = {OuterMargin("0pt", "0pt","40pt","40pt")};
                                Tabla32.Width = "500pt";
                                Tabla32.Border="solid";
                                Tabla32.ColSep="solid";
                                Tabla32.RowSep="solid";
                                Tabla32.HAlign="center";
                                Tabla32.TableEntriesHAlign="center";
                                Tabla32.Style= [Tabla32.Style {NumberFormat("%.2f")}];
                                add(s32,Tabla32);
                    
                                add(chap3,s3);
                                add(chap3,s32);
                                add(rpt,chap3);
                                
                                %%%%%%%Capitulo5
                    
                    
                                chap5 = Chapter('Factores de seguridad estaticos');

                                 if R==1
                                     s5=Section("Factor de seguridad al desplazamiento con suelo posterior");
                                     image5 = mlreportgen.report.FormalImage();
                                     image5.Image = which('fsdespo.png');
                                     image5.Caption = 'Ecucaion para el Factor de seguridad al deslizamiento con suelo posterior ';
                                     add(s5,image5);                                 

                                 end

                                 if R==0

                                     s5=Section("Factor de seguridad al desplazamiento");
                                     image5 = mlreportgen.report.FormalImage();
                                     image5.Image = which('fsdes.png');
                                     image5.Caption = 'Ecucaion para el Factor de seguridad al deslizamiento ';
                                     add(s5,image5);
    %                     
    
    %                                
                        

                                 end
                                s52=Section("Factor de seguridad al volcamiento");
%                                 plot52 = Image("fsvol.png");
%                                 plot52.Style = {OuterMargin("100pt", "0pt","0pt","5pt")};
%                                 append(s52,plot52);
                                image52 = mlreportgen.report.FormalImage();
                                image52.Image = which('fsvol.png');
                                image5.Caption = 'Ecucaion para el Factor de seguridad al volcamiento ';
                                add(s52,image52);
                                
                                s53=Section("Resultados");
                                %TABLARESULTADOS
                                fsde=[FSD]';
                                fsv1=[FSVR]';
                                fsv2=[FSVA]';
                                
                               
                                Tabla5= FormalTable(["F.S DES","F.S.1 VOL","F.S.2 VOL"],[fsde,fsv1,fsv2]);
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