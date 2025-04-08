#03 script for semantic plots

library(ggrepel)

#condition_half_thaw_vs_control.down-over   45    -none- character
#condition_full_freeze_vs_control.up-over  357    -none- character
#condition_full_thaw_vs_control.up-over    220    -none- character
#condition_half_freeze_vs_control.up-over    7    -none- character (no filter)
#condition_full_freeze_vs_control.up-under  34    -none- character (no filter)

# half-thaw up-over, 45 total terms ----
#note the revgio website converts the adjusted pvaleus to log and stores them in the column "value". Use this to plot the size of the circles

# Combine the data for the three categories into one data frame
revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");

#BP
revigo.data_1 <- rbind(c("GO:0002262","myeloid cell homeostasis",0.029223643500396353,4.974169607335022,-3.3680055973462504,4.07459704459004,-7.238824186844268,0.5231576076878256,0.87798065),
                       c("GO:0002376","immune system process",0.9082220067499159,-2.7154588844376977,-3.0088590239328323,5.567019304402755,-4.555955204081924,1,-0),
                       c("GO:0002520","immune system development",0.03717135215556184,5.389744520187939,0.6372133130834228,4.17906322239498,-7.931814138253839,0.6630459491888464,0.54582388),
                       c("GO:0006811","monoatomic ion transport",4.917982551075166,-1.4844699800329935,-6.475856054425679,6.300613307421874,-4.665546248849069,0.9684763948745858,0.26008441),
                       c("GO:0015669","gas transport",0.06168396612174118,-2.3105472648646095,-6.1233865547708,4.3990157256487645,-4.089375595110798,0.9730258725552657,0.26422016),
                       c("GO:0015701","bicarbonate transport",0.05227913652391295,-0.684651229018001,-6.932227837517204,4.327174958937619,-5.696803942579511,0.9732731359446819,-0),
                       c("GO:0030097","hemopoiesis",0.1392485814361512,5.825427595811932,3.037201222601514,4.7526245626267665,-7.407823242604133,0.7339219987121813,0.57215337),
                       c("GO:0030099","myeloid cell differentiation",0.05429006187890528,5.3914320017525545,2.9138283411691948,4.343566132386124,-6.614393726401688,0.6970468032131346,0.78692007),
                       c("GO:0030221","basophil differentiation",0.00022890582376289568,5.213175812502073,3.575275687539176,1.9731278535996986,-4.049148541111453,0.7482771315433124,0.6781901),
                       c("GO:0033280","response to vitamin D",0.0014817344720996044,-5.458833603476449,-3.0527677066762324,2.780317312140151,-4.014573525916998,0.9737546155888971,0.43592814),
                       c("GO:0042592","homeostatic process",1.9296957851447603,1.8862235130926577,7.973027929563176,5.894315508737042,-5.889410289700751,1,-0),
                       c("GO:0044058","regulation of digestive system process",0.007706496066684156,-3.196387459216515,4.757272350873843,3.495821753385906,-4.031984286006359,0.8788938400644237,0.62753179),
                       c("GO:0044706","multi-multicellular organism process",0.02605342090892743,6.479083317115236,-2.5246261480480325,4.024731889655249,-4.113509274827518,0.8237870659067879,0.36836569),
                       c("GO:0045597","positive regulation of cell differentiation",0.1878332271915899,-4.089770102810914,3.0096607065600014,4.882604217705935,-5.2373214362725635,0.8533351140685633,0.50390783),
                       c("GO:0045639","positive regulation of myeloid cell differentiation",0.018179552842072556,-3.86990911535464,3.834947837622317,3.8684680990209865,-4.591760034688151,0.8205738831457021,0.76700703),
                       c("GO:0046697","decidualization",0.0018287852371594788,6.173272309719724,0.7062458401200971,2.8715729355458786,-4.518557371497695,0.7398630238885812,0.32959467),
                       c("GO:0048513","animal organ development",1.042075302533505,5.9854446230745495,2.660636589552154,5.626726235442988,-4.906578314837765,0.7477811707514167,0.59413566),
                       c("GO:0048534","hematopoietic or lymphoid organ development",0.014758272250347556,5.7103605222407845,0.4749988090369488,3.7779340488377793,-8.276544327964814,0.6721172692571121,0.55731491),
                       c("GO:0048821","erythrocyte development",0.010379525363528291,5.213493613825126,-0.6180233878038863,3.6251065754034677,-9.3829996588791,0.40134867311214034,0),
                       c("GO:0048872","homeostasis of number of cells",0.04988670146264914,4.949247784257959,-3.5990957254845237,4.306832322684375,-6.291579099865287,0.5382118913608908,0.81223771),
                       c("GO:0048878","chemical homeostasis",1.4930591569935086,4.072166956866794,-5.532745723507018,5.782903837869686,-6.015022873584507,0.5959087698175691,0.87809953),
                       c("GO:0050878","regulation of body fluid levels",0.10579633357850049,-1.1449713080648396,6.532646379261836,4.633306827560638,-4.798602875679548,0.9509804120526695,0.17736222),
                       c("GO:0051240","positive regulation of multicellular organismal process",0.38622565851159724,-3.4230893093622883,3.964094936874705,5.195669996526704,-5.425968732272281,0.8353878706987813,0.2027298),
                       c("GO:0061515","myeloid cell development",0.015257926897916028,5.1959124576673155,3.1010671920553357,3.792391689498254,-8.455931955649724,0.6991485993591323,0.83676157),
                       c("GO:0065008","regulation of biological quality",3.028571729559731,-4.724035542840641,1.5589714159574155,6.0900643235623235,-7.169411331314856,0.9391820316512361,-0),
                       c("GO:0098771","inorganic ion homeostasis",0.9913542024229459,3.8458175043415275,-5.50252484399147,5.605056036739983,-7.723538195826756,0.5969020726517068,0.56163674),
                       c("GO:1903708","positive regulation of hemopoiesis",0.03034602044271765,-3.680291846376105,3.665024120551252,4.090963076595732,-4.399027104313252,0.8162823218391655,0.88248958),
                       c("GO:1904374","cellular response to kainic acid",2.461352943687051E-06,-5.645052615360363,-2.568012703663678,0.3010299956639812,-4.120330794367947,0.9717134168654659,0.00560173))

#CC
revigo.data_2 <- rbind(c("GO:0005833","hemoglobin complex",0.012910797096992051,-3.3933698640412597,5.560651863176492,3.7139103541289553,-5.156144577376839,0.9263863700589318,0.14995624),
                       c("GO:0005886","plasma membrane",16.678688713171073,5.625754852283576,1.0790619572995892,6.8250353366409255,-5.4089353929735005,0.856619111938595,0.29772069),
                       c("GO:0031838","haptoglobin-hemoglobin complex",0.012267003967687075,-5.237662989010902,0.8061614091117291,3.6917002082901615,-6.137272471682025,0.9264273094419068,0),
                       c("GO:0032589","neuron projection membrane",0.05668373924531724,5.16769087181656,3.454382138438226,4.356350978030569,-4.1732774798310075,0.8566525118503909,-0),
                       c("GO:0071944","cell periphery",18.431522806558757,0.018401395326788084,-5.725350811367457,6.868434645240301,-5.258848401148215,0.9999158775841537,3.744E-05))


#MF
revigo.data_3 <- rbind(c("GO:0004900","erythropoietin receptor activity",0.0008512683013680994,5.352046285027899,0.04691755237571196,2.5921767573958667,-4.120330794367947,0.9864477495024467,0.00822618),
                       c("GO:0005200","structural constituent of cytoskeleton",0.14988652047857953,0.6017447026058276,6.3853436574071765,4.836767047394205,-6.6595558851598815,0.9818543765362813,0),
                       c("GO:0015267","channel activity",1.673927439592528,-5.844372943130174,-0.762892504925089,5.88473533996407,-4.13430394008393,0.7508376461095972,0.01411561),
                       c("GO:0019215","intermediate filament binding",0.008097962559168332,0.7781268886138406,-6.419959008261827,3.569490954348783,-3.9802068210508446,1,-0),
                       c("GO:0022803","passive transmembrane transporter activity",1.673927439592528,-5.86290473851858,0.5480040989232389,5.88473533996407,-4.114638779968488,0.7508376461095972,0.51036347))

# Combine both data frames and label them with categories
revigo.data_1 <- data.frame(revigo.data_1)
revigo.data_2 <- data.frame(revigo.data_2)
revigo.data_3 <- data.frame(revigo.data_3)

names(revigo.data_1) <- revigo.names
names(revigo.data_2) <- revigo.names
names(revigo.data_3) <- revigo.names
revigo.data_1$Category <- "Biological process"
revigo.data_2$Category <- "Cellular component"
revigo.data_3$Category <- "Molecular function"

# Combine both data frames
combined_data <- rbind(revigo.data_1, revigo.data_2, revigo.data_3)

# Clean the data (convert to numeric)
combined_data$plot_X <- as.numeric(as.character(combined_data$plot_X))
combined_data$plot_Y <- as.numeric(as.character(combined_data$plot_Y))
combined_data$log_size <- as.numeric(as.character(combined_data$log_size))
combined_data$value <- as.numeric(as.character(combined_data$value))
combined_data$frequency <- as.numeric(as.character(combined_data$frequency))
combined_data$uniqueness <- as.numeric(as.character(combined_data$uniqueness))
combined_data$dispensability <- as.numeric(as.character(combined_data$dispensability))

# Make three faceted plots. Note we are only displaying the "cluster represantative" term 

# Filter data for specific descriptions (cluster representative)
filtered_data <- subset(combined_data, description %in% c(
  "immune system process", 
  "bicarbonate transport", 
  "homeostatic process", 
  "erythrocyte development", 
  "regulation of biological quality", 
  "cellular response to kainic acid", 
  "cell periphery", 
  "neuron projection membrane", 
  "haptoglobin-hemoglobin complex", 
  "intermediate filament binding", 
  "channel activity", 
  "structural constituent of cytoskeleton", 
  "erythropoietin receptor activity"
))

# Create faceted plots with reduced size legend breaks
p_faceted <- ggplot(data = combined_data) +
  geom_point(aes(plot_X, plot_Y, colour = Category, size = value), alpha = 0.6) +
  scale_colour_manual(values = c("Cellular component" = "#fbb4ae", 
                                 "Biological process" = "#b3cde3", 
                                 "Molecular function" = "#ccebc5")) +
  # Reduce number of size legend breaks
  scale_size_continuous(range = c(5, 15), breaks = c(-9.3, -6.3, -4)) +  
  geom_point(aes(plot_X, plot_Y, size = value), shape = 21, 
             fill = "transparent", colour = I(adjustcolor("black", alpha.f = 0.6))) +
  # Add labels only for filtered data
  geom_label_repel(data = filtered_data, aes(plot_X, plot_Y, label = description), 
                   fill = "white", colour = I(adjustcolor("black", alpha.f = 0.85)), 
                   size = 5, nudge_y = -1, 
                   max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3) +
  labs(y = "semantic space x", x = "semantic space y", size = "log(adjusted p value)") +
  theme_bw() +
  theme(legend.key = element_blank(),
        legend.text = element_text(size = 18),  # Increase the size of legend text
        legend.title = element_text(size = 20), # Increase the size of legend titles
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        strip.text = element_text(size = 22),
        plot.margin = margin(10, 10, 10, 10),         # Add margins around the plot
        panel.spacing = unit(1.5, "lines"),
        axis.text.x = element_text(margin = margin(t = 5))  # Add extra space at the top of x-axis labels
  ) +
  xlim(min(combined_data$plot_X) - 1, max(combined_data$plot_X) + 1) +  # Add buffer
  ylim(min(combined_data$plot_Y) - 1, max(combined_data$plot_Y) + 1) +  # Add buffer +
  facet_wrap(~Category, scales = "free") +  # Separate plots by Category
  guides(colour = guide_legend(override.aes = list(size = 16, shape = 21, fill = c("#fbb4ae", "#b3cde3", "#ccebc5"), colour = "black")))  # Add black outline to legend

p_faceted

ggsave("semantic_half_thaw.pdf", p_faceted, width = 32, height = 8, units = "in" )


# half-freeze up-over, 7 total terms ----
#note the revgio website converts the adjusted pvaleus to log and stores them in the column "value". Use this to plot the size of the circles

# Combine the data for the two categories into one data frame
revigo_data_1 <- rbind(
  c("GO:0005771", "multivesicular body", 0.08128512088026983, -0.4271008517240548, -0.19648531995175533, 4.512897756320646, -5.865198684182015, 0.2720362940678587, 0),
  c("GO:0031906", "late endosome lumen", 5.240176633877717E-05, -1.9866091505364973, 0.1260169747266498, 1.3424226808222062, -5.786706221391004, 0.2631255761129791, 0.52418255),
  c("GO:0045334", "clathrin-coated endocytic vesicle", 0.0304354449540031, 0.7262897921587944, -1.6782436714403381, 4.086288629021682, -5.77758113552173, 0.38975251631010027, 0.48084031),
  c("GO:0097208", "alveolar lamellar body", 0.002188397575195599, 1.5553601138424256, 0.8284290555583542, 2.9434945159061026, -5.118141188227927, 0.4332869891616213, 0.40990714)
)
revigo_data_2 <- rbind(
  c("GO:0036296", "response to increased oxygen levels", 0.002431816708362806, 0, 1.5453809621887546, 2.9951962915971793, -5.390082710976574, 1, 0),
  c("GO:0055093", "response to hyperoxia", 0.0017180243546935615, 1.5453809621887546, 0, 2.8444771757456815, -5.5221582696250495, 1, 0.64708963)
)

# Combine both data frames and label them with categories
revigo_data_1 <- data.frame(revigo_data_1)
revigo_data_2 <- data.frame(revigo_data_2)
names(revigo_data_1) <- revigo.names
names(revigo_data_2) <- revigo.names
revigo_data_1$Category <- "Cellular component"
revigo_data_2$Category <- "Biological process"

# Combine both data frames
combined_data <- rbind(revigo_data_1, revigo_data_2)

# Clean the data (convert to numeric)
combined_data$plot_X <- as.numeric(as.character(combined_data$plot_X))
combined_data$plot_Y <- as.numeric(as.character(combined_data$plot_Y))
combined_data$log_size <- as.numeric(as.character(combined_data$log_size))
combined_data$value <- as.numeric(as.character(combined_data$value))
combined_data$frequency <- as.numeric(as.character(combined_data$frequency))
combined_data$uniqueness <- as.numeric(as.character(combined_data$uniqueness))
combined_data$dispensability <- as.numeric(as.character(combined_data$dispensability))

p_combined <- ggplot(data = combined_data) +
  geom_point(aes(plot_X, plot_Y, colour = Category, size = value), alpha = 0.6) +
  scale_colour_manual(values = c("Cellular component" = "#fbb4ae", 
                                 "Biological process" = "#b3cde3" 
  )) +
  # Reduce number of size legend breaks
  scale_size_continuous(range = c(5, 15)) +  
  geom_point(aes(plot_X, plot_Y, size = value), shape = 21, 
             fill = "transparent", colour = I(adjustcolor("black", alpha.f = 0.6))) +
  # Add labels for filtered data
  geom_label(data = combined_data, aes(plot_X, plot_Y, label = description), 
             fill = "white", colour = I(adjustcolor("black", alpha.f = 0.85)), 
             size = 4, nudge_y = -0.3) +
  labs(y = "semantic space x", x = "semantic space y", size = "log(adjusted p value)") +
  theme_bw() +
  theme(legend.key = element_blank(),
        legend.text = element_text(size = 18),  # Increase the size of legend text
        legend.title = element_text(size = 20), # Increase the size of legend titles
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        strip.text = element_text(size = 22),
        plot.margin = margin(10, 10, 10, 10),         # Add margins around the plot
        panel.spacing = unit(1.5, "lines"),
        axis.text.x = element_text(margin = margin(t = 5))) +  # Add extra space at the top of x-axis labels
  xlim(min(combined_data$plot_X) - 1, max(combined_data$plot_X) + 1) +  # Add buffer to x-axis
  ylim(min(combined_data$plot_Y) - 1, max(combined_data$plot_Y) + 1) +  # Add buffer to y-axis
  facet_wrap(~Category, scales = "free") +  # Separate plots by Category with independent axes
  guides(colour = guide_legend(override.aes = list(size = 16, shape = 21, 
                                                   fill = c("#fbb4ae", "#b3cde3"), 
                                                   colour = "black")))  # Add black outline to legend

# Print the plot
print(p_combined)

ggsave("semantic_half_freeze.pdf", p_combined, width = 16, height = 6, units = "in" )



# full-freeze, up-under, 34 total ----
#Molecular function

revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data_1 <- rbind(c("GO:0003676","nucleic acid binding",20.48386614107179,-2.037333887571077,-1.0234750511899808,6.972410120534942,-7.787812395596042,0.6529160259876793,0),
                     c("GO:0003723","RNA binding",6.057692897451659,-5.4255169775949454,-0.5255899180296434,6.443305554056935,-5.621602099051862,0.7054309078186742,0.26095929),
                     c("GO:0097159","organic cyclic compound binding",39.06621935354991,2.3956414250921063,-4.205169549961636,7.252799553663059,-6.282329496997738,0.711883544075481,0.17598293),
                     c("GO:1901363","heterocyclic compound binding",19.934158544093723,4.57940596824329,1.1514370823525253,6.960596100738509,-6.006123085058789,0.7432177856612713,0.1389612));

#Cellular component
revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data_2 <- rbind(c("GO:0005654","nucleoplasm",1.7423412635088946,-4.5836721208374005,-3.959991515356944,5.844007212503244,-7.481486060122113,0.4100458874736729,0.81200181),
                     c("GO:0005694","chromosome",2.4935006213976125,-4.841724412128465,1.7787179701804874,5.99968328386563,-4.423658649794207,0.732613293688829,0.18986071),
                     c("GO:0031974","membrane-enclosed lumen",4.185416413755364,1.968176738743003,5.6798070923166035,6.224612300524052,-5.787812395596042,0.9999079530169853,5.19E-05),
                     c("GO:0031981","nuclear lumen",3.2078514814166494,-4.922970425555403,-3.3315115448725288,6.109087959211095,-6.903089986991944,0.37596347624410836,0),
                     c("GO:0043227","membrane-bounded organelle",30.260550315864123,-5.668977592740207,-0.5818644012650819,7.083750226085466,-5.268411234813262,0.614842953911049,0.54229438),
                     c("GO:0043231","intracellular membrane-bounded organelle",29.492520084536526,-4.869466955596324,-1.2689098004391817,7.072585289060935,-5.853871964321762,0.5880726019243736,0.3577277),
                     c("GO:1902494","catalytic complex",6.912374390158849,4.8817356940823755,-2.6185966711149136,6.442500777747844,-6.928117992693875,0.8994442536544398,0.46448686),
                     c("GO:1990904","ribonucleoprotein complex",4.3863098139649965,4.489355533778273,-4.026232235005877,6.2449729190612935,-5.357535479757878,0.8994442536544398,-0));


#Biological process
revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data_3 <- rbind(c("GO:0006259","DNA metabolic process",5.572970721566783,-3.3047694962943734,5.927254101818165,6.3549130597720795,-4.124360062995832,0.7176167665548577,0.64503033),
                     c("GO:0006325","chromatin organization",1.6774440287109929,0.2707573693389707,-5.557897509249895,5.833474781755866,-4.567030709125595,0.9362052124779853,0.44379736),
                     c("GO:0006396","RNA processing",4.402491558667012,-3.8338562585197327,5.887423192039941,6.2525248812610785,-4.856985199745905,0.6904292040899399,0.02047605),
                     c("GO:0008152","metabolic process",58.08643542977757,4.336316008237624,-4.638952369733295,7.372900851027942,-4.756961951313706,1,-0),
                     c("GO:0009889","regulation of biosynthetic process",14.06877837293839,5.1230562122704715,3.7974614494277947,6.757082570711665,-4.510041520575165,0.6424638985624395,0.79032866),
                     c("GO:0009892","negative regulation of metabolic process",2.9251604469836634,4.918250092243592,4.461345231210766,6.0749761643750535,-4.137868620686963,0.7160717183489771,-0),
                     c("GO:0010467","gene expression",13.02581206063927,-4.673249810924311,3.7263998442728594,6.723630996379145,-6.440093374963888,0.7566250959934206,0.58931936),
                     c("GO:0010468","regulation of gene expression",13.833180130521608,4.994191546773767,3.310776106846075,6.7497482160021365,-4.47755576649368,0.6279776326404726,0.56504465),
                     c("GO:0010556","regulation of macromolecule biosynthetic process",13.904945798300691,5.411771072937269,3.389556465530805,6.751995483723061,-4.478861916295964,0.6276953551023128,0.85807699),
                     c("GO:0016070","RNA metabolic process",8.105171248384922,-3.8382002807574627,5.171660080833773,6.517588433595825,-4.3655227298392685,0.7019043792166101,0.6180467),
                     c("GO:0031326","regulation of cellular biosynthetic process",13.947256455402673,5.282583596680055,3.006295714696119,6.7533149697461,-4.083019952679618,0.6280233474933207,0.85374817),
                     c("GO:0033044","regulation of chromosome organization",0.1861028960721779,5.525905814288523,0.8118363817900471,4.878584981900468,-4.078313524516398,0.8563803047377999,0.19734623),
                     c("GO:0043170","macromolecule metabolic process",33.33895178165279,-5.162614881567509,-1.931412509229938,7.131778079538736,-8.931814138253838,0.8654431770859236,0.10787326),
                     c("GO:0044237","cellular metabolic process",39.640441131550894,-4.268010458168969,-0.9293586012882235,7.206964612305403,-7.863279432843593,0.8516904820688888,0.27679073),
                     c("GO:0044238","primary metabolic process",48.680252087830326,-5.477605179056666,-0.1020780580517672,7.2961789470216925,-8.067526235322847,0.8504626014192497,0.2504376),
                     c("GO:0051276","chromosome organization",1.530548023678806,-0.4033205193277124,-5.658783380212555,5.793673765853863,-6.790484985457369,0.9363449667749678,0),
                     c("GO:0090304","nucleic acid metabolic process",13.388603177774023,-4.032329462666487,4.384235476187825,6.7355614561517365,-7.617982957425132,0.7073375150731304,0.53830305));


# Combine both data frames and label them with categories
revigo.data_1 <- data.frame(revigo.data_1)
revigo.data_2 <- data.frame(revigo.data_2)
revigo.data_3 <- data.frame(revigo.data_3)

names(revigo.data_1) <- revigo.names
names(revigo.data_2) <- revigo.names
names(revigo.data_3) <- revigo.names
revigo.data_1$Category <- "Molecular function"
revigo.data_2$Category <- "Cellular component"
revigo.data_3$Category <- "Biological process"

# Combine both data frames
combined_data <- rbind(revigo.data_1, revigo.data_2, revigo.data_3)

# Clean the data (convert to numeric)
combined_data$plot_X <- as.numeric(as.character(combined_data$plot_X))
combined_data$plot_Y <- as.numeric(as.character(combined_data$plot_Y))
combined_data$log_size <- as.numeric(as.character(combined_data$log_size))
combined_data$value <- as.numeric(as.character(combined_data$value))
combined_data$frequency <- as.numeric(as.character(combined_data$frequency))
combined_data$uniqueness <- as.numeric(as.character(combined_data$uniqueness))
combined_data$dispensability <- as.numeric(as.character(combined_data$dispensability))

# Make three faceted plots. Note we are only displaying the "cluster represantative" term 

# Filter data for specific descriptions (cluster representative)
filtered_data <- subset(combined_data, description %in% c(
  "nucleic acid binding", "RNA binding", "organic cyclic compound binding", "heterocyclic compound binding", "membrane-enclosed lumen", "nuclear lumen", "ribonucleoprotein complex", "RNA processing", "metabolic process", "negative regulation of metabolic process", "chromosome organization", "regulation of chromosome organization", "primary metabolic process"
))

# Create faceted plots with reduced size legend breaks
p_faceted <- ggplot(data = combined_data) +
  geom_point(aes(plot_X, plot_Y, colour = Category, size = value), alpha = 0.6) +
  scale_colour_manual(values = c("Cellular component" = "#fbb4ae", 
                                 "Biological process" = "#b3cde3", 
                                 "Molecular function" = "#ccebc5")) +
  # Reduce number of size legend breaks
  scale_size_continuous(range = c(5, 15), breaks = c(-8.9, -7.5, -4.1)) +  
  geom_point(aes(plot_X, plot_Y, size = value), shape = 21, 
             fill = "transparent", colour = I(adjustcolor("black", alpha.f = 0.6))) +
  # Add labels only for filtered data
  geom_label_repel(data = filtered_data, aes(plot_X, plot_Y, label = description), 
                   fill = "white", colour = I(adjustcolor("black", alpha.f = 0.85)), 
                   size = 5, nudge_y = -1, 
                   max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3) +
  labs(y = "semantic space x", x = "semantic space y", size = "log(adjusted p value)") +
  theme_bw() +
  theme(legend.key = element_blank(),
        legend.text = element_text(size = 18),  # Increase the size of legend text
        legend.title = element_text(size = 20), # Increase the size of legend titles
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        strip.text = element_text(size = 22),
        plot.margin = margin(10, 10, 10, 10),         # Add margins around the plot
        panel.spacing = unit(1.5, "lines"),
        axis.text.x = element_text(margin = margin(t = 5))) + # Add extra space at the top of x-axis labels
  xlim(min(combined_data$plot_X) - 1, max(combined_data$plot_X) + 1) +  # Add buffer
  ylim(min(combined_data$plot_Y) - 1, max(combined_data$plot_Y) + 1) +  # Add buffer +
  facet_wrap(~Category, scales = "free") +  # Separate plots by Category
  guides(colour = guide_legend(override.aes = list(size = 16, shape = 21, fill = c("#fbb4ae", "#b3cde3", "#ccebc5"), colour = "black")))  # Add black outline to legend

p_faceted

ggsave("semantic_full_freeze_up_under.pdf", p_faceted, width = 32, height = 8, units = "in")


# full-thaw, up-over, 220 ---- 

#BP
revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data_1 <- rbind(c("GO:0000122","negative regulation of transcription by RNA polymerase II",0.48144309713813077,-5.351990312465718,-1.0096907667922088,5.291373291067754,-5.078313524516398,0.6260230188777682,0.66006709),
                     c("GO:0001501","skeletal system development",0.12859092318998627,1.9297698211389995,7.474275508016173,4.718044733513918,-3.522430206021097,0.9716035202497083,0.62260371),
                     c("GO:0001704","formation of primary germ layer",0.027832979087213168,2.101274165599062,7.372129421264694,4.053424204069344,-3.728562684472197,0.974205641431024,0.38313698),
                     c("GO:0001944","vasculature development",0.14776978532719576,2.016787801314285,7.422186539535101,4.778418982771813,-3.570727484753387,0.9713694885524307,0.63131261),
                     c("GO:0002573","myeloid leukocyte differentiation",0.025118106790326354,2.3264319116804772,7.231481443005882,4.008855563996213,-3.364589647111321,0.9726303156987943,0.43744914),
                     c("GO:0002684","positive regulation of immune system process",0.33138425357330603,-3.4456644656833393,2.1276576645491314,5.129161200385081,-3.430863788251155,0.669520832084355,0.81261334),
                     c("GO:0002931","response to ischemia",0.0029806984148050184,4.716105263367756,-3.5805526069709384,3.0835026198302673,-6.104577453960592,0.9005615491052023,0.16229824),
                     c("GO:0006351","DNA-templated transcription",2.0416135034942102,2.145637857932235,-7.322084386054139,5.918800159642457,-4.1090204030103115,0.9729638271284617,0.55749194),
                     c("GO:0006355","regulation of DNA-templated transcription",11.135305937063894,-4.832368368124002,0.5147569920852789,6.655528356088988,-3.735130253580378,0.6328357732658609,0.72663885),
                     c("GO:0006357","regulation of transcription by RNA polymerase II",4.1876991326832185,-5.086343023335698,0.6061934023883713,6.230801833816229,-4.987162775294828,0.6697415474267465,0.52820302),
                     c("GO:0006366","transcription by RNA polymerase II",0.5182403236462522,2.290889986989241,-7.245177781223761,5.323359371070169,-4.924453038607469,0.9761978230033163,0.30819511),
                     c("GO:0006950","response to stress",6.596819705552286,5.213317637623641,-0.7064962727068468,6.428160883324095,-4.590066876668706,0.8448269799121558,0.42938362),
                     c("GO:0006952","defense response",1.1245700077941203,5.211612021528977,-1.3528850795293175,5.6598135537959715,-3.824569370762764,0.8482915931580096,0.54266122),
                     c("GO:0006955","immune response",0.6332667307635791,5.655002916146227,-1.1119322227987707,5.410414467099877,-3.852728794025814,0.8762227112180321,0.32014509),
                     c("GO:0006979","response to oxidative stress",0.7667902052527142,5.281921658169942,-1.4916433091402408,5.493504057283511,-3.9036114533126334,0.8532569346265914,0.33211207),
                     c("GO:0007339","binding of sperm to zona pellucida",0.019688362196552716,0.45788747442771427,-8.119922116935903,3.9030899869919438,-5.51427857351842,0.9682013333242276,0.01012068),
                     c("GO:0008285","negative regulation of cell population proliferation",0.12504657495107693,-5.590076142558383,-2.5104918273316486,4.705906455700312,-3.5774971590641984,0.7125317329832288,0.49141741),
                     c("GO:0009266","response to temperature stimulus",0.3946804058731623,5.266028813254126,-1.9509838523891527,5.205074381062382,-3.729568071643146,0.8505735307780258,0.82180635),
                     c("GO:0009314","response to radiation",0.322284631740495,5.2621196971680835,-2.051765788457453,5.117069019829141,-5.521433504406157,0.8528856923436969,0.6868663),
                     c("GO:0009605","response to external stimulus",1.5774441813148754,5.018165068391862,-0.8529026970638248,5.806780784629541,-7.080921907623926,0.8656009381257714,0.24783982),
                     c("GO:0009607","response to biotic stimulus",0.9712178739906423,5.600308047963087,-0.8114692544285643,5.596143873584628,-6.198596289982645,0.8714621468372578,0.33573219),
                     c("GO:0009612","response to mechanical stimulus",0.04310075139690395,4.957220006209572,-2.5256679916994726,4.2433357485598755,-9.962573502059376,0.8525994913932721,0),
                     c("GO:0009617","response to bacterium",0.20962358480205134,4.707954647973894,-1.6995649774081967,4.930271349530451,-5.653647025549361,0.8617016171890152,0.85305276),
                     c("GO:0009628","response to abiotic stimulus",0.9321586641272724,5.378583556243717,-0.9759536778446853,5.578317093840565,-5.354577730650908,0.8719349135281987,0.33417055),
                     c("GO:0009636","response to toxic substance",1.0264457113410923,5.340922206187146,0.25728885217156905,5.620163132442645,-5.1830961606243395,0.8171143618072587,0.46373282),
                     c("GO:0009725","response to hormone",0.6238151354598208,5.2127538899540395,0.5660793060070904,5.403883727818924,-4.133712660915805,0.8090608346422051,0.62354679),
                     c("GO:0009889","regulation of biosynthetic process",14.06877837293839,-4.865420495945428,0.38227712554833276,6.757082570711665,-3.465879933957705,0.6536003884445853,0.74907494),
                     c("GO:0009892","negative regulation of metabolic process",2.9251604469836634,-5.087562867482097,-0.8881581044144531,6.0749761643750535,-4.09799710864927,0.6067544184891552,0.66866099),
                     c("GO:0009968","negative regulation of signal transduction",0.6342832695293219,-4.213014198128402,-2.4551354390043056,5.411111047994698,-4.187755303199631,0.5862917053211523,0.69943835),
                     c("GO:0009988","cell-cell recognition",0.0283597086171622,0.531738148416301,-8.092241486127268,4.061565561884839,-4.705533773838407,0.9796745033139246,0.81092229),
                     c("GO:0010286","heat acclimation",0.007019778595395469,4.993711819663435,-2.9171697429251817,3.4553017716570764,-5.749579997691106,0.8713282773672214,0.58756768),
                     c("GO:0010466","negative regulation of peptidase activity",0.1552375301583423,-6.191559648233925,-0.9128209414913679,4.7998297168509225,-4.140861702705469,0.6128632063257,0.84111908),
                     c("GO:0010556","regulation of macromolecule biosynthetic process",13.904945798300691,-4.661175561216705,0.406329402632083,6.751995483723061,-3.48286218787106,0.6363352463415168,0.81341281),
                     c("GO:0010558","negative regulation of macromolecule biosynthetic process",2.375589561717219,-5.0001428473450895,-0.7737154565591675,5.984597964833603,-3.894738541320597,0.5846673246244126,0.86710949),
                     c("GO:0010628","positive regulation of gene expression",0.484610858376656,-4.56646617884702,1.7586822557359176,5.294221453198817,-4.818156412055227,0.6438752668375276,0.75098073),
                     c("GO:0010648","negative regulation of cell communication",0.6498513268981425,-4.428806489035791,-2.203489317802612,5.421641761483437,-4.540607512240769,0.6233052649773564,0.76421433),
                     c("GO:0010803","regulation of tumor necrosis factor-mediated signaling pathway",0.008407981655634965,-3.0147370385704524,-4.814224840282429,3.5336449787987627,-5.341035157335565,0.7291532394036275,0.8357282),
                     c("GO:0010821","regulation of mitochondrion organization",0.04492707528111973,-2.4168415680958364,4.62258025508094,4.2613580461941245,-3.556643121181756,0.7446422511834158,0.55120798),
                     c("GO:0010968","regulation of microtubule nucleation",0.005518353299746368,-2.3550433048642097,5.700480567983122,3.350829273582968,-6.522878745280337,0.7357199785669883,0.86394122),
                     c("GO:0014074","response to purine-containing compound",0.022218633022663005,5.3011079116566515,1.4625909522845515,3.9555915504057246,-5.404503778174425,0.8400239816407576,0.70302803),
                     c("GO:0016070","RNA metabolic process",8.105171248384922,1.9809718644911136,-7.408828219037334,6.517588433595825,-4.045757490560675,0.9798407668174134,0.44797897),
                     c("GO:0019219","regulation of nucleobase-containing compound metabolic process",11.820748427527752,-5.010649869839059,0.4304751973765862,6.681471171615664,-3.468623281918542,0.6534786473840323,0.76741413),
                     c("GO:0019222","regulation of metabolic process",15.757534779778567,-3.8908295614624353,0.2809832192288051,6.806314448734751,-3.71263950451249,0.7201189719264951,0.27459697),
                     c("GO:0023057","negative regulation of signaling",0.6509343221933648,-4.557775700095272,-2.512989628966347,5.422364920152503,-4.525783735923745,0.6306463062327422,0.76366091),
                     c("GO:0030512","negative regulation of transforming growth factor beta receptor signaling pathway",0.013217465307599461,-4.40209047164883,-3.824097560286892,3.7300551523755,-3.9609906176452157,0.6739525709615822,0.71519724),
                     c("GO:0031110","regulation of microtubule polymerization or depolymerization",0.07398088542840169,-2.516609388745804,4.356680038579672,4.477960080114487,-3.47514191490446,0.7146343185642694,0.8660324),
                     c("GO:0031112","positive regulation of microtubule polymerization or depolymerization",0.01683073142893205,-3.3549399186046367,4.163781058808398,3.83499260373303,-4.774690718274138,0.668267361606496,0.76733883),
                     c("GO:0031113","regulation of microtubule polymerization",0.018521680901245056,-2.523096910491397,5.109300429927535,3.8765642139838454,-4.2020403562628035,0.7205455580404194,0.83359284),
                     c("GO:0031326","regulation of cellular biosynthetic process",13.947256455402673,-4.59937412864384,0.4230183506218002,6.7533149697461,-3.638538487024328,0.6383073876199139,0.85486863),
                     c("GO:0031334","positive regulation of protein-containing complex assembly",0.06646145218543774,-3.4737970761215218,3.546587641481522,4.431412016420789,-5.253365801062421,0.6536569781443842,0.72731537),
                     c("GO:0031397","negative regulation of protein ubiquitination",0.021199632903976567,-6.509713170117028,-2.032020479081475,3.9352048674265814,-3.838039365694781,0.6946466830208851,0.69187457),
                     c("GO:0031400","negative regulation of protein modification process",0.14815129503346727,-5.951112048230418,-1.4812403200946835,4.779538773870285,-4.690369832574101,0.6569598670438831,0.6874985),
                     c("GO:0031960","response to corticosteroid",0.01687503578191842,5.290808061699787,1.5594421481791392,3.8361341494653747,-4.467245621007502,0.8339388069257687,0.83174918),
                     c("GO:0032069","regulation of nuclease activity",0.011260689717368256,-7.77574876060594,1.2792884768552018,3.6604860157849677,-5.43533393574791,0.757227331073158,0.59101894),
                     c("GO:0032075","positive regulation of nuclease activity",0.00044796623575104325,-6.4410224160156115,2.873084621503585,2.2624510897304293,-6.486782399932061,0.7205474125632312,0.83999138),
                     c("GO:0032273","positive regulation of protein polymerization",0.04537258016392709,-3.4538340916619017,3.7287367083358003,4.265643141942135,-3.8535680645369492,0.6589860239274896,0.84300737),
                     c("GO:0032677","regulation of interleukin-8 production",0.009426981774321403,-6.316247268623095,3.2474933087769706,3.5833121519830775,-5.785156151952302,0.7727490089020477,0.69623331),
                     c("GO:0032757","positive regulation of interleukin-8 production",0.007667114419585163,-4.980388850968603,3.0606734089752923,3.493597449000527,-6.2168113089247425,0.7054925099188839,0.52133071),
                     c("GO:0032774","RNA biosynthetic process",6.5238873564778945,2.199627376687657,-7.2932774254777435,6.423332724148457,-4.052566278112949,0.9708310554448416,0.73359902),
                     c("GO:0033120","positive regulation of RNA splicing",0.037476559920579036,-4.962776694273487,2.4522118570925304,4.18261434773635,-6.151195298948196,0.6991285602295557,0.61084849),
                     c("GO:0033993","response to lipid",0.26565382321214337,5.420798123904634,0.4765491110360334,5.0331462008952395,-3.418036695510142,0.8360715618960046,0.57289911),
                     c("GO:0034599","cellular response to oxidative stress",0.33220388410355384,5.259809277865052,-0.02879712872300819,5.1302340301621,-4.329754146925876,0.8062030693266209,0.58535014),
                     c("GO:0034605","cellular response to heat",0.1474424253856854,5.14263736656517,-2.169001664917309,4.7774558227219925,-4.950781977329818,0.8380775398778346,0.75578826),
                     c("GO:0034620","cellular response to unfolded protein",0.11508547958797542,4.871088451079562,-2.077544354670973,4.669855926622827,-3.7106808932462694,0.8566935105443331,0.43637235),
                     c("GO:0035914","skeletal muscle cell differentiation",0.011287764599748814,2.4417910621974768,7.157670263253454,3.6615287401319825,-5.278189384787454,0.9738332349520927,0.00975343),
                     c("GO:0042026","protein refolding",0.19468809513975832,2.5653028173101324,-7.102152530495428,4.8981709930144,-6.432973633840939,0.9847245744310757,-0),
                     c("GO:0042127","regulation of cell population proliferation",0.401650957409684,-2.2553021697835383,0.46209756432243104,5.212677574440172,-4.931814138253839,0.794502095534774,0.18065212),
                     c("GO:0043086","negative regulation of catalytic activity",0.32922810839463623,-7.402660528635388,1.0032292455313516,5.126326260094941,-3.6344427581495715,0.7488505616094537,0.89709315),
                     c("GO:0043154","negative regulation of cysteine-type endopeptidase activity involved in apoptotic process",0.01250121160098653,-6.110213240595173,-1.8005803817472976,3.7058637122839193,-5.785156151952302,0.6282612132593653,0.8455894),
                     c("GO:0043254","regulation of protein-containing complex assembly",0.27547216010451103,-2.684908418651534,3.396013767472652,5.048907701483771,-3.504797127144315,0.7029916595690094,0.84247536),
                     c("GO:0043487","regulation of RNA stability",0.25213361149247043,-5.840423716479934,0.6493367224401623,5.010461090711631,-3.7585705755479037,0.7219869142846689,0.79208545),
                     c("GO:0043488","regulation of mRNA stability",0.23649663624122655,-5.907563246803413,0.6035804131021125,4.982655594477174,-3.945409924595048,0.7222603833727552,0.29525168),
                     c("GO:0044089","positive regulation of cellular component biogenesis",0.17646916065058676,-3.7689352941230063,2.956245086195429,4.855500983970402,-3.573662816625531,0.6691830582368015,0.75475725),
                     c("GO:0044092","negative regulation of molecular function",0.41366974383370786,-7.260525065147119,1.1598154986466331,5.225482447973504,-3.6965836102354577,0.7542878401054559,0.59400974),
                     c("GO:0045089","positive regulation of innate immune response",0.0929874528595531,-2.9519398332336984,-1.0603951168379508,4.577261953585815,-3.4270761462592803,0.6091743276004301,0.83813902),
                     c("GO:0045471","response to ethanol",0.004935012652092536,5.048554599315185,1.9756659224220021,3.302330928684399,-3.8050030685818315,0.8667200950183707,0.48431185),
                     c("GO:0045597","positive regulation of cell differentiation",0.1878332271915899,-3.5910109164529627,2.4396179608332575,4.882604217705935,-3.765136127875084,0.6841753861621416,0.80936936),
                     c("GO:0045639","positive regulation of myeloid cell differentiation",0.018179552842072556,-4.280019893829771,3.5460461711595257,3.8684680990209865,-6.55129368009492,0.6949088434519853,0.50711865),
                     c("GO:0045646","regulation of erythrocyte differentiation",0.004405821769199821,-0.6873961896763594,2.9330435765661282,3.2530955858490316,-5.841637507904751,0.7710869979944814,0.85253065),
                     c("GO:0045648","positive regulation of erythrocyte differentiation",0.0026533384732946408,-4.513558321087924,4.512731689926842,3.0330214446829107,-5.147520006363144,0.7248556331376669,0.82995461),
                     c("GO:0045934","negative regulation of nucleobase-containing compound metabolic process",1.5483509895204945,-5.305048429371428,-0.8861670163797098,5.798696212904897,-3.9896873200754137,0.6042798516209292,0.75901503),
                     c("GO:0045935","positive regulation of nucleobase-containing compound metabolic process",1.8080311091383094,-4.6146606007301605,1.480327301481553,5.86603259646637,-5.146910470148135,0.6119869276816542,0.85423234),
                     c("GO:0046677","response to antibiotic",0.32127301568063965,5.136962287826737,0.4764845969312681,5.115703683637,-4.458420756053419,0.8336430803696918,0.58345237),
                     c("GO:0046683","response to organophosphorus",0.022378620964002664,4.994821146562044,1.2545422363391145,3.9587071910872846,-5.698970004336019,0.8622525668414808,0.33918578),
                     c("GO:0048545","response to steroid hormone",0.0706703657191426,5.294903008593819,1.2174978930856584,4.458078570949283,-4.0400051616715835,0.8217217308056273,0.79689335),
                     c("GO:0048585","negative regulation of response to stimulus",0.7695395364908127,-4.285884010901243,-2.618847531249022,5.495058433016688,-4.318758762624412,0.6107087650187871,0.7445616),
                     c("GO:0048856","anatomical structure development",3.2670128346264704,2.2004307954995532,7.30980035330893,6.122977274352911,-3.3137812472578236,0.9687836850894861,0.54051877),
                     c("GO:0050896","response to stimulus",16.94347616187058,3.685982414357038,6.35636346477149,6.83782868559364,-4.667561540084395,1,-0),
                     c("GO:0051051","negative regulation of transport",0.1163653831186927,-5.754554216648765,-2.475041573615992,4.67465909631365,-3.529785088486167,0.6955905297071991,0.60939339),
                     c("GO:0051091","positive regulation of DNA-binding transcription factor activity",0.0783867071976015,-6.729873754199271,1.0463313667942324,4.503082164576042,-3.535422529502235,0.7215200707084825,0.65994658),
                     c("GO:0051094","positive regulation of developmental process",0.3250364243315372,-3.7960072358759067,2.480238046430105,5.12076142698027,-3.4438459085353963,0.6945584228593574,0.60968307),
                     c("GO:0051131","chaperone-mediated protein complex assembly",0.03827157692138995,-0.5719868007607728,-8.338220990790584,4.191730393362857,-6.158015195409886,0.9905489310409717,0.01059744),
                     c("GO:0051248","negative regulation of protein metabolic process",0.5470258463226723,-5.570651920577084,-1.181598335311544,5.34683590735996,-5.206209615309182,0.634652183186493,0.14621276),
                     c("GO:0051254","positive regulation of RNA metabolic process",1.7391009199503533,-4.679741505133983,1.5177757334494901,5.849151506011876,-5.876148359032914,0.6089118276187235,0.42178087),
                     c("GO:0051384","response to glucocorticoid",0.015917569486824157,5.213256480441686,1.5932561361481634,3.8107700112343634,-4.732828271596986,0.8316643610774099,0.72766113),
                     c("GO:0051385","response to mineralocorticoid",0.0013512827660841907,5.206481861549779,2.1394250743972556,2.7403626894942437,-3.358132426439526,0.8529746060949714,0.87041758),
                     c("GO:0051412","response to corticosterone",0.0005168841181742806,5.128513061211802,2.2872151597732793,2.3242824552976926,-4.024568191490737,0.8541831456688163,0.83043054),
                     c("GO:0051495","positive regulation of cytoskeleton organization",0.08146093702426663,-3.4340599263540583,3.4085251429143217,4.519788629954198,-3.595303090589227,0.6504635830586895,0.78538336),
                     c("GO:0051591","response to cAMP",0.008725496185370593,5.213009174123134,1.747229399290175,3.549738731264899,-6.430626090384954,0.8448871815918273,0.17249493),
                     c("GO:0051592","response to calcium ion",0.07172874748492802,5.4093622758608575,0.8851171805047412,4.464534256300944,-3.8451775231773766,0.8426604498451692,0.81489482),
                     c("GO:0051707","response to other organism",0.9426071073732241,4.891692976429256,-1.3244359232221965,5.58315795063656,-6.357535479757878,0.8439967540322622,0.65797034),
                     c("GO:0051716","cellular response to stimulus",13.811589142499587,5.251297353303597,-0.618846115260681,6.749069834706852,-4.252588192113577,0.8294472601897692,0.6276778),
                     c("GO:0060700","regulation of ribonuclease activity",0.008885484126710252,-7.780114086536703,1.165322953866152,3.5576274884268266,-7.329754146925876,0.7437557408269576,0.67717689),
                     c("GO:0060759","regulation of response to cytokine stimulus",0.023286860200223185,-3.192283529484266,-5.028353939739139,3.9759829437125465,-3.3665608350359806,0.7529289076912541,0.35663899),
                     c("GO:0070370","cellular heat acclimation",0.005638959593987033,5.058832411872164,-2.849703464627952,3.3602146132953523,-8.30539480106643,0.8688367084306275,0.52836855),
                     c("GO:0070424","regulation of nucleotide-binding domain, leucine rich repeat containing receptor signaling pathway",0.0020380002373728777,-2.5698488006285904,-4.474131393726434,2.9185545305502734,-7.178486471595227,0.712041915627155,0.73109203),
                     c("GO:0070426","positive regulation of nucleotide-binding domain, leucine rich repeat containing receptor signaling pathway",0.0005414976476111512,-2.3176480548563467,-2.5514378239291498,2.3443922736851106,-8.707743928643524,0.6755734257083764,0.88215161),
                     c("GO:0070434","positive regulation of nucleotide-binding oligomerization domain containing 2 signaling pathway",0.00042827541220154676,-2.30002065155178,-2.7025306961534667,2.2430380486862944,-8.707743928643524,0.6789661063103913,0.19956512),
                     c("GO:0070887","cellular response to chemical stimulus",2.1752846591599666,5.26499589193062,0.12068081302141091,5.946342695005116,-3.3787686556689507,0.8011546147761546,0.7161887),
                     c("GO:0071248","cellular response to metal ion",0.09907437868929116,5.378633138449635,0.7107935872437287,4.604798253272669,-3.3450816506715206,0.8340003802749788,0.84262182),
                     c("GO:0071277","cellular response to calcium ion",0.06120154094477852,5.472295750358438,0.7990705445739695,4.395605929313231,-4.151810883008602,0.8393501543250926,0.50268364),
                     c("GO:0071496","cellular response to external stimulus",0.008789491361906458,4.462297890335182,-3.178288236322747,3.552911450216509,-3.8408570826546673,0.8911040397214242,0.47796468),
                     c("GO:0090063","positive regulation of microtubule nucleation",0.004164609180718489,-3.252161143307065,4.691232528773687,3.2286569581089353,-7.184422251675732,0.6813848233926025,0.30040949),
                     c("GO:0090083","regulation of inclusion body assembly",0.002660722532125702,-2.223358232664453,6.104092398028347,3.0342272607705505,-6.229884705212898,0.7749878755714237,0.49984055),
                     c("GO:0090084","negative regulation of inclusion body assembly",0.0008196305302477877,-6.252635667134342,-2.426727658664861,2.5237464668115646,-7.127843727251707,0.7174332777238509,0.47283998),
                     c("GO:0090224","regulation of spindle organization",0.02457414778977151,-2.4106104958742307,5.01896914586592,3.9993480692067216,-4.0254883072626715,0.7328773871990751,0.79061027),
                     c("GO:0090559","regulation of membrane permeability",0.020360311550179283,-0.5719323378409427,-5.874844389978142,3.917663024327375,-3.600383380750477,0.8342854899318428,0.46625749),
                     c("GO:0097305","response to alcohol",0.06957260230625817,5.130854172395576,1.0024281238390274,4.451279718904047,-5.0467236633326955,0.845621550564132,0.5614304),
                     c("GO:1900034","regulation of cellular response to heat",0.0020724591785844966,-2.935405549353243,-5.610075672938915,2.9258275746247424,-4.399027104313252,0.7607119846339153,0.66333501),
                     c("GO:1901029","negative regulation of mitochondrial outer membrane permeabilization involved in apoptotic signaling pathway",0.008218457478971061,-3.6107015019423647,-2.202776107294845,3.5237464668115646,-6.467245621007502,0.558455509049284,0.51050828),
                     c("GO:1901099","negative regulation of signal transduction in absence of ligand",0.009643580833365865,-4.456889973090432,-3.9924449077059116,3.5931752634781025,-5.809668301829708,0.6824933221907967,0.53913954),
                     c("GO:1901673","regulation of mitotic spindle assembly",0.009591892421548436,-2.4056922973596113,5.444428974973914,3.5908418347816027,-5.628932137728264,0.7329299052749255,0.69505169),
                     c("GO:1901700","response to oxygen-containing compound",0.9029891703916372,5.22910296320948,0.2881404360844113,5.564509832144328,-3.980360378055731,0.8190935128873925,0.64838158),
                     c("GO:1902235","regulation of endoplasmic reticulum stress-induced intrinsic apoptotic signaling pathway",0.006357674653543651,-3.3783165839583944,-4.82362102030813,3.4122925093230463,-6.17134010346468,0.7039591452436795,0.71155363),
                     c("GO:1902236","negative regulation of endoplasmic reticulum stress-induced intrinsic apoptotic signaling pathway",0.002872398885282788,-4.187494734457519,-4.42209448816739,3.0674428427763805,-5.690369832574102,0.6642341386363274,0.82340891),
                     c("GO:1902380","positive regulation of endoribonuclease activity",2.461352943687051E-06,-6.525685823881169,3.8558004333819063,0.3010299956639812,-8.707743928643524,0.7750921733578381,-0),
                     c("GO:1902905","positive regulation of supramolecular fiber organization",0.07551430831231871,-3.5129387586651375,3.4602720798310127,4.486869510668216,-3.871675058870552,0.6616323067937543,0.7588537),
                     c("GO:1903265","positive regulation of tumor necrosis factor-mediated signaling pathway",0.0022029108845999104,-2.68354398610423,-2.6707693681352795,2.9523080096621253,-6.863279432843593,0.6841412085831138,0.44459572),
                     c("GO:1903533","regulation of protein targeting",0.026484157674072667,-0.6375040659995526,4.401338608680195,4.031852631395629,-4.313363730737707,0.7939180868570836,0.71744557),
                     c("GO:1903573","negative regulation of response to endoplasmic reticulum stress",0.00882148895017439,-4.213941544589669,-4.086790821526933,3.554489160003819,-4.543633966870957,0.6710216270329745,0.87203781),
                     c("GO:1903706","regulation of hemopoiesis",0.06241744929895991,-1.9833707868960877,1.4867945319158613,4.404149249209695,-5.1605219526258015,0.7410447896941746,0.75435881),
                     c("GO:1903708","positive regulation of hemopoiesis",0.03034602044271765,-4.114798043927851,3.197703558664582,4.090963076595732,-5.767003889607846,0.6849364484913475,0.88248958),
                     c("GO:1903747","regulation of establishment of protein localization to mitochondrion",0.0026311862968014573,-0.5055021574644858,6.4350147178128205,3.0293837776852097,-4.804100347590766,0.818891645174938,0.62262435),
                     c("GO:1903845","negative regulation of cellular response to transforming growth factor beta stimulus",0.00021659905904446047,-4.1984602119520975,-5.44875874769339,1.9493900066449128,-3.9418132070888463,0.7455291871018658,0.75578202),
                     c("GO:1903955","positive regulation of protein targeting to mitochondrion",0.0015161934133112232,-4.247556485251409,5.294313015885219,2.7902851640332416,-7.4698003017969175,0.748153806999686,0.27127975),
                     c("GO:1905897","regulation of response to endoplasmic reticulum stress",0.02520179279041171,-2.9796177916387716,-4.325172821134265,4.010299956639812,-4.497572880015567,0.7264029048692053,0.51970642),
                     c("GO:2000116","regulation of cysteine-type endopeptidase activity",0.05986502629635645,-7.200031895238542,0.7048329325591243,4.386017139806847,-3.390321301008853,0.7103169195084599,0.87486517),
                     c("GO:2000117","negative regulation of cysteine-type endopeptidase activity",0.019272393549069605,-6.709244269521325,-1.2550411318929895,3.893817223967463,-5.515700160653214,0.6570541456932089,0.52409253),
                     c("GO:2001141","regulation of RNA biosynthetic process",11.138259560596317,-4.7796182441793995,0.49967293437943994,6.655643536761422,-3.5984926221309,0.6336430632910381,0.86932043),
                     c("GO:2001234","negative regulation of apoptotic signaling pathway",0.06547198830207555,-4.240358362369886,-3.2542258417858196,4.4248979631844,-3.6414103954917336,0.6215906145818052,0.88962256),
                     c("GO:2001236","regulation of extrinsic apoptotic signaling pathway",0.03807713003883867,-3.4161245287042252,-3.962721275460306,4.189518386126378,-3.4927015157145633,0.6946906431685244,0.89243889),
                     c("GO:2001239","regulation of extrinsic apoptotic signaling pathway in absence of ligand",0.01037214130469723,-3.425691855430236,-4.61667218395455,3.624797578960761,-4.910094888560602,0.71519477965044,0.85451985),
                     c("GO:2001240","negative regulation of extrinsic apoptotic signaling pathway in absence of ligand",0.009643580833365865,-4.2459946058921645,-3.946051080737863,3.5931752634781025,-5.809668301829708,0.6583803678141219,0.78555021),
                     c("GO:2001242","regulation of intrinsic apoptotic signaling pathway",0.03628772644877819,-3.3437215338239814,-3.9591753626160826,4.168615322211018,-4.217527375833714,0.6925290614774259,0.79037321),
                     c("GO:2001243","negative regulation of intrinsic apoptotic signaling pathway",0.021536838257261694,-4.218760593898737,-3.6224943364504054,3.942057683841395,-4.534617148551582,0.6401645723830887,0.87466986));

#CC
revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data_2 <- rbind(c("GO:0002199","zona pellucida receptor complex",0.0027373684606494552,-2.9773603355769196,5.39335555913722,3.040602340114073,-6.3979400086720375,0.9531418007871624,0.11855093),
                     c("GO:0016234","inclusion body",0.015009363072749745,-0.8963969215588212,-6.297620082266439,3.779307827583586,-3.6020391456931224,0.9564612310904073,0.028543),
                     c("GO:0016235","aggresome",0.01093450190935817,-5.605791769294568,-1.1576081211892322,3.6417714706539592,-4.987162775294828,0.9570021056673127,0.02398244),
                     c("GO:0016607","nuclear speck",0.1457268168659327,5.382773389705159,0.4064993523094479,4.7664202835980785,-3.5797652697763676,0.795050699860043,0.15920637),
                     c("GO:0035976","transcription factor AP-1 complex",0.0004965691191150789,3.050478090327195,5.297835832158311,2.3010299956639813,-8.838631997765026,0.8481537758078661,0),
                     c("GO:1904813","ficolin-1-rich granule lumen",0.0003094199536194461,5.761548547631324,-2.2470529819658425,2.0969100130080562,-3.5233494859996863,0.8343416112642571,0.31646148))

#MF
revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data_3 <- rbind(c("GO:0000976","transcription cis-regulatory region binding",3.1067997063691872,4.551468522240395,3.722483490699583,6.153311705691067,-4.847711655616943,0.5910246034011546,0.56045648),
                     c("GO:0000978","RNA polymerase II cis-regulatory region sequence-specific DNA binding",1.5700290520409343,4.229215502166939,3.926039775942141,5.856906437642981,-4.51427857351842,0.6167235095413515,0.8821991),
                     c("GO:0000981","DNA-binding transcription factor activity, RNA polymerase II-specific",2.318396477687505,0.17486838381781214,-5.73738377654272,6.026186262324628,-4.886056647693163,0.7980411056860798,0.59502457),
                     c("GO:0000987","cis-regulatory region sequence-specific DNA binding",1.6055225747295157,4.543208408643803,4.085142895271948,5.866615151983271,-4.421360790031928,0.6159174015766371,0.88479294),
                     c("GO:0001067","transcription regulatory region nucleic acid binding",3.1069437671586497,5.235858602477053,2.3509372821671692,6.153331843233675,-4.209714835966758,0.776390441201567,0.38192782),
                     c("GO:0001228","DNA-binding transcription activator activity, RNA polymerase II-specific",0.1573820470090874,0.6479769854424573,-5.498623183600779,4.857959358058422,-8.367542707815275,0.8348652420391767,-0),
                     c("GO:0002020","protease binding",0.038785093453871176,-5.5761164306249285,2.718778687424243,4.249687427805301,-3.528335777129688,0.8398425614096409,0.3240228),
                     c("GO:0003690","double-stranded DNA binding",3.5854285829222525,5.096311427859945,3.3331392826616932,6.215539485865457,-4.0726296369609765,0.6978154779970837,0.57288986),
                     c("GO:0003700","DNA-binding transcription factor activity",5.948455529733023,-0.039043417527101595,-5.37871493992275,6.435402523670657,-4.359518563029578,0.7923286097926779,0.76677956),
                     c("GO:0003712","transcription coregulator activity",0.684238546943507,-0.4477611560137711,-5.7578969475856105,5.496207067291576,-4.856985199745905,0.7974106747782107,0.44150046),
                     c("GO:0003714","transcription corepressor activity",0.2346401022063281,-0.7234471434982792,-5.560597383869107,5.0314044242841645,-4.97469413473523,0.8113461907542548,0.79583041),
                     c("GO:0003796","lysozyme activity",0.03902082929117311,6.861711916733293,-2.4913743135701942,4.2523189329416535,-3.4133586931543953,0.997564849157584,0.0280213),
                     c("GO:0008134","transcription factor binding",0.19164887024800478,-4.466945967257158,2.2238931619113957,4.94350935487179,-4.430626090384954,0.8254649986094608,0.23989505),
                     c("GO:0031072","heat shock protein binding",0.16150960690110552,-4.222533297300894,2.8990105123838483,4.869202374517823,-3.809982228048801,0.8270892664947646,0.36045995),
                     c("GO:0031249","denatured protein binding",9.385778707391866E-05,-3.7025551857121353,4.748578558551349,1.6434526764861874,-8.707743928643524,0.8792317518061382,0),
                     c("GO:0043565","sequence-specific DNA binding",4.369488160527723,4.548088544984891,3.066290335487484,6.301428929011041,-5.2564902352715706,0.6913461804476809,0.02612876),
                     c("GO:0044183","protein folding chaperone",0.39308732596174273,-5.995916086658217,-4.464054929864955,5.2554895980755365,-5.576754126063192,1,-0),
                     c("GO:0051059","NF-kappaB binding",0.013768719089820438,-4.9904732218724615,3.7095367966437878,3.7999605274059833,-4.826813731587726,0.8481837040652798,0.30186974),
                     c("GO:0055131","C3HC4-type RING finger domain binding",0.0002750251435189244,-5.641935326847255,0.7089617609042863,2.103803720955957,-7.396855627379818,0.84375954900473,0.17137154),
                     c("GO:0061783","peptidoglycan muralytic activity",0.1852949163370239,-0.08074051859139668,8.081677962102708,4.928866765396401,-3.2760184031292554,0.997564849157584,-0),
                     c("GO:0097718","disordered domain specific binding",0.007506440226679215,-5.401565812330029,1.2882446598379222,3.53655844257153,-4.889410289700751,0.8232324760018942,0.57555242),
                     c("GO:0140110","transcription regulator activity",6.680930431014821,5.491906231921531,-4.5514474553280495,6.48583523646154,-6.29929628285498,1,-0))


# Combine both data frames and label them with categories
revigo.data_1 <- data.frame(revigo.data_1)
revigo.data_2 <- data.frame(revigo.data_2)
revigo.data_3 <- data.frame(revigo.data_3)

names(revigo.data_1) <- revigo.names
names(revigo.data_2) <- revigo.names
names(revigo.data_3) <- revigo.names
revigo.data_1$Category <- "Biological process"
revigo.data_2$Category <- "Cellular component"
revigo.data_3$Category <- "Molecular function"

# Combine both data frames
combined_data <- rbind(revigo.data_1, revigo.data_2, revigo.data_3)

# Clean the data (convert to numeric)
combined_data$plot_X <- as.numeric(as.character(combined_data$plot_X))
combined_data$plot_Y <- as.numeric(as.character(combined_data$plot_Y))
combined_data$log_size <- as.numeric(as.character(combined_data$log_size))
combined_data$value <- as.numeric(as.character(combined_data$value))
combined_data$frequency <- as.numeric(as.character(combined_data$frequency))
combined_data$uniqueness <- as.numeric(as.character(combined_data$uniqueness))
combined_data$dispensability <- as.numeric(as.character(combined_data$dispensability))

# Make three faceted plots. Note we are only displaying the "cluster represantative" term 

# Filter data for specific descriptions (cluster representative)
filtered_data <- subset(combined_data, description %in% c("nucleic acid binding", "RNA binding", "organic cyclic compound binding", "heterocyclic compound binding", "membrane-enclosed lumen", "nuclear lumen", "ribonucleoprotein complex", "RNA processing", "metabolic process", "negative regulation of metabolic process", "chromosome organization", "regulation of chromosome organization", "primary metabolic process", "response to ischemia", "response to stress", "immune response", "response to oxidative stress", "response to temperature stimulus", "heat acclimation", "regulation of metabolic process", "negative regulation of protein ubiquitination", "response to corticosteroid", "skeletal muscle cell differentiation", "protein refolding", "response to stimulus", "chaperone-mediated protein complex assembly", "positive regulation of endoribonuclease activity", "chromosome organization", "nucleic acid binding", "lysozyme activity", "DNA-binding transcription activator activity", "RNA polymerase II-specific", "heat shock protein binding", "denatured protein binding", "sequence-specific DNA binding", "protein folding chaperone", "peptidoglycan muralytic activity", "transcription regulator activity", "NF-kappaB binding", "zona pellucida receptor complex", "inclusion body", "aggresome", "nuclear speck", "transcription factor AP-1 complex", "ficolin-1-rich granule lumen", "DNA-binding transcription activator activity, RNA polymerase II-specific", "response to heat"))

# Ensure 'Category' is a factor with levels ordered as desired
combined_data$Category <- factor(combined_data$Category, 
                                 levels = c("Biological process", "Cellular component", "Molecular function"))

# Create faceted plots with adjusted layout
p_faceted <- ggplot(data = combined_data) +
  geom_point(aes(plot_X, plot_Y, colour = Category, size = value), alpha = 0.6) +
  scale_colour_manual(values = c("Cellular component" = "#fbb4ae", 
                                 "Biological process" = "#b3cde3", 
                                 "Molecular function" = "#ccebc5")) +
  # Reduce number of size legend breaks
  scale_size_continuous(range = c(5, 15), breaks = c(-9.9, -6.5, -3.3)) +  
  geom_point(aes(plot_X, plot_Y, size = value), shape = 21, 
             fill = "transparent", colour = I(adjustcolor("black", alpha.f = 0.6))) +
  # Add labels only for filtered data
  geom_label_repel(data = filtered_data, aes(plot_X, plot_Y, label = description), 
                   fill = "white", colour = I(adjustcolor("black", alpha.f = 0.85)), 
                   size = 5, nudge_y = -1, 
                   max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3) +
  labs(y = "semantic space x", x = "semantic space y", size = "log(adjusted p value)") +
  theme_bw() +
  theme(legend.key = element_blank(),
        legend.text = element_text(size = 18),  # Increase the size of legend text
        legend.title = element_text(size = 20), # Increase the size of legend titles
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        strip.text = element_text(size = 22),
        plot.margin = margin(10, 10, 10, 10),         # Add margins around the plot
        panel.spacing = unit(1.5, "lines"),
        axis.text.x = element_text(margin = margin(t = 5))  # Add extra space at the top of x-axis labels
  ) +
  xlim(min(combined_data$plot_X) - 1, max(combined_data$plot_X) + 1) +  # Add buffer
  ylim(min(combined_data$plot_Y) - 1, max(combined_data$plot_Y) + 1) +  # Add buffer
  facet_wrap(~Category, scales = "free", nrow = 1) +  # Two rows for facets
  guides(colour = guide_legend(override.aes = list(size = 16, shape = 21, fill = c("#fbb4ae", "#b3cde3", "#ccebc5"), colour = "black")))  # Add black outline to legend

p_faceted
# Save the plot
ggsave("semantic_full_thaw.pdf", p_faceted, width = 32, height = 8, units = "in")




# full-freeze, up-over, 357 ----

revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data_1 <- rbind(c("GO:0001508","action potential",0.06162489365109269,-4.799161010878112,-5.428489034403624,4.398599635095352,-6.790484985457369,0.8220805078447191,0.68478641),
                     c("GO:0003008","system process",1.1200115821424117,0.46838510335087497,-7.874275825765757,5.658049574713654,-15.543633966870956,0.8187018703672141,0.5744028),
                     c("GO:0003012","muscle system process",0.12514995177471178,0.19949308791404607,-7.911136231491089,4.706265334409592,-29.80966830182971,0.8157325221140677,0),
                     c("GO:0003013","circulatory system process",0.15395024256879394,0.33014236505618244,-7.983743414692848,4.796213427493061,-13.320572103387882,0.8133925087968915,0.6750183),
                     c("GO:0006090","pyruvate metabolic process",0.6101743174459072,0.6740988139473534,6.10900087044685,5.3942817824175275,-5.034798298974088,0.9243104981000739,0.01507009),
                     c("GO:0006164","purine nucleotide biosynthetic process",1.420257259625133,1.2767029495048938,5.647163181411971,5.761193877021937,-3.9894710812450622,0.8630309085050462,0.61726269),
                     c("GO:0006936","muscle contraction",0.11758129147287409,0.2145868360368148,-8.047684291845746,4.679173423631723,-24.546681659952963,0.7831175237826169,0.66147192),
                     c("GO:0007010","cytoskeleton organization",1.8978237258769566,5.60453519783376,-1.560401719494347,5.8870825415010435,-9.442492798094342,0.9512750699117933,0.47793149),
                     c("GO:0007154","cell communication",9.236864411598072,5.800191133628669,-4.797303210931867,6.574350790896292,-3.3791178497146355,0.9908173491369955,0.02075352),
                     c("GO:0007165","signal transduction",8.851839893322994,-4.5848671541248684,-0.068499182556585,6.555859777363872,-3.754734160542539,0.7833224517486792,0.48667797),
                     c("GO:0007271","synaptic transmission, cholinergic",0.03497582532979299,-1.6861152794211098,-0.4488303021694152,4.152624639447619,-3.0925004992900154,0.8732924271252744,0.30022974),
                     c("GO:0007386","compartment pattern specification",0.00046765705930053964,1.5813850358223247,-7.153255222090761,2.2810333672477277,-3.163493594842893,0.8116055486660203,0.35846492),
                     c("GO:0008016","regulation of heart contraction",0.05387409323142216,-7.854875146851097,-1.4760595085734358,4.34022592125719,-17.732828271596986,0.7638194875595071,-0),
                     c("GO:0009187","cyclic nucleotide metabolic process",0.160918332752372,1.4493376687904314,5.735875688111058,4.815438273573759,-3.1886998661854187,0.9008240218781088,0.28242909),
                     c("GO:0009611","response to wounding",0.15748228404298487,-3.2384131805712837,6.350611848298152,4.806064599188238,-3.223649007228483,0.9769583046655569,0.22790045),
                     c("GO:0009612","response to mechanical stimulus",0.04310075139690395,-2.6862345647238746,6.359093036883728,4.2433357485598755,-5.023191662661934,0.9724175355167151,0.20715612),
                     c("GO:0009653","anatomical structure morphogenesis",1.4683422507330033,2.830045571943146,-6.118030813013707,5.775654130082595,-7.777283528852417,0.7910527372416032,0.66427229),
                     c("GO:0010522","regulation of calcium ion transport into cytosol",0.0009993092951369425,-3.964781888239433,-6.174762015326029,2.60959440922522,-4.22767829327708,0.7379068849095035,0.69105487),
                     c("GO:0010649","regulation of cell communication by electrical coupling",0.0006916401771760612,-5.157688226603912,4.696327821043564,2.450249108319361,-6.1487416512809245,0.8894380906749374,0.0998374),
                     c("GO:0010650","positive regulation of cell communication by electrical coupling",1.7229470605809354E-05,-5.6019012595088,5.18755952224053,0.9030899869919435,-3.2585537168714382,0.8803726880579987,0.2810339),
                     c("GO:0010752","regulation of cGMP-mediated signaling",0.008550740126368814,-5.328474991396612,3.5921723673784225,3.540954808926133,-3.174398641792348,0.8708750863683778,0.50149267),
                     c("GO:0010817","regulation of hormone levels",0.28198736134645064,-4.996556487617959,-4.705754159543851,5.059059541001648,-3.178111947695968,0.8231458176890047,0.59351112),
                     c("GO:0014706","striated muscle tissue development",0.08190890326001768,3.1152540767969,-6.083847218058563,4.522170267708623,-16.906578314837766,0.794425537201223,0.49769722),
                     c("GO:0014728","regulation of the force of skeletal muscle contraction",0.0013315919425346944,-8.865386803896351,-1.218282098858195,2.733999286538387,-3.247867007502326,0.8020738358649028,0.6920104),
                     c("GO:0016525","negative regulation of angiogenesis",0.02570390879092387,-7.032242832650592,0.6815594686351238,4.018866863150907,-3.610421088992206,0.7553615381296211,0.64601059),
                     c("GO:0019722","calcium-mediated signaling",0.1524365105084264,-3.250438834023171,1.0449241604224064,4.791922117501735,-5.692503962086787,0.8389608938144636,0.17216733),
                     c("GO:0022603","regulation of anatomical structure morphogenesis",0.9122389347540132,-6.413710151581105,0.021546761566275757,5.568935882496253,-5.24184537803261,0.8030327862225702,0.2061066),
                     c("GO:0023052","signaling",9.184767415191992,-5.154162407031687,-1.1887392103147802,6.571894386105714,-4.411168274405793,0.8241509868737548,0.31189867),
                     c("GO:0030029","actin filament-based process",0.8043086081733359,5.340322202729704,1.322206021586638,5.514250152357573,-18.237321436272563,0.992772260410035,-0),
                     c("GO:0030048","actin filament-based movement",0.08138217373006863,5.08395915270661,-3.2063586657547067,4.51936852704,-11.288192770958808,0.9228583202203658,0.69308715),
                     c("GO:0031032","actomyosin structure organization",0.12413341300896903,5.326345497096468,-2.2734306061486143,4.702723414104734,-15.638272163982407,0.9093433802309,0.01298598),
                     c("GO:0031033","myosin filament organization",0.005744797770565576,5.32943167624522,-2.4535728577040103,3.368286884902131,-8.536107011014092,0.9241896118366025,0.58575866),
                     c("GO:0031034","myosin filament assembly",0.0016564905311013851,5.101362719543666,-2.2704831060973505,2.82865989653532,-8.536107011014092,0.9232058816665253,0.54460305),
                     c("GO:0032409","regulation of transporter activity",0.15078248133026873,-7.323409241012887,-3.5390250587756467,4.787184081777398,-6.511449283499556,0.8517174803229269,0.1412233),
                     c("GO:0032501","multicellular organismal process",3.5834862120669206,5.92983985485175,2.201102516095032,6.163132142633867,-4.175223537524454,1,-0),
                     c("GO:0032502","developmental process",3.830665120083749,1.975936052455447,2.711415023121784,6.192100572841233,-6.06398920428479,1,-0),
                     c("GO:0032835","glomerulus development",0.012124624600602412,1.643925796732653,-6.907227562845972,3.692582562274909,-3.267281445046938,0.7602954036694877,0.52427255),
                     c("GO:0032879","regulation of localization",0.9058074195121588,-6.025250619138435,-1.9402709661876212,5.565863160327171,-8.48412615628832,0.8574042808314767,0.18646309),
                     c("GO:0032970","regulation of actin filament-based process",0.3388347689338468,-3.4399624974732506,-1.0569814405094489,5.138817229617316,-8.818156412055227,0.8631649177698836,0.17197811),
                     c("GO:0032989","cellular anatomical entity morphogenesis",0.07956569525762759,5.73394054378892,-0.9983222934754977,4.509565403222792,-7.552841968657781,0.9683356085075253,0.26533795),
                     c("GO:0035051","cardiocyte differentiation",0.027151184321811853,1.9100813750679946,-6.837328316368108,4.042654253167793,-10.892790030352131,0.7145887730127507,0.45954329),
                     c("GO:0035150","regulation of tube size",0.07517710295903358,-5.114852761176296,-5.28067686966679,4.484925911049592,-4.257274868695302,0.8346859326075463,0.53196899),
                     c("GO:0035994","response to muscle stretch",0.0030250027677913853,-2.6619630383796884,6.560950382691505,3.089905111439398,-4.826813731587726,0.9711445369698556,0.51024383),
                     c("GO:0036309","protein localization to M-band",1.9690823549496406E-05,2.7892550053869893,0.024065445923598775,0.9542425094393249,-3.251546049054529,0.9946200139584319,0.11950981),
                     c("GO:0040011","locomotion",0.4992706765092561,3.3869829484721183,3.267937204052964,5.307164307090027,-4.079876673709276,1,-0),
                     c("GO:0040012","regulation of locomotion",0.4079643277102412,-6.52467709612531,-2.5609405827477523,5.219450916891037,-6.040481623027001,0.866027740024571,0.18985527),
                     c("GO:0042391","regulation of membrane potential",0.3881799727488847,-5.12510278769669,-4.506290398260443,5.197861985482968,-3.8639266014706206,0.8191405566182615,0.61058821),
                     c("GO:0042592","homeostatic process",1.9296957851447603,6.338198650885205,0.9617391011254881,5.894315508737042,-3.559212096705457,1,-0),
                     c("GO:0042692","muscle cell differentiation",0.13121718678090036,2.9678496616915764,-5.850510754499195,4.726824975390675,-24.282329496997736,0.7757003983089082,0.51598038),
                     c("GO:0042866","pyruvate biosynthetic process",0.003148070414975738,0.3481606059073544,6.2823294331276776,3.1072099696478683,-3.5751804098697746,0.9388317956850951,0.56590066),
                     c("GO:0044057","regulation of system process",0.15059049580066114,-7.467533272884338,-1.3942652494784464,4.786630768031185,-14.244887733604928,0.7810336572203173,0.67239884),
                     c("GO:0045663","positive regulation of myoblast differentiation",0.009284223303587556,-7.195177934732642,1.8706122629195843,3.5766868052009957,-3.468699946924814,0.8154639661849677,0.5989596),
                     c("GO:0045844","positive regulation of striated muscle tissue development",0.007467744831146512,-7.9015644538684615,0.5315712940874417,3.482158695411276,-3.453251342262446,0.7652458189500077,0.6079081),
                     c("GO:0046068","cGMP metabolic process",0.056625885822464285,1.4779060872406293,5.507621310567841,4.361859992489294,-3.0654687254138993,0.8952100972781292,0.62088878),
                     c("GO:0046365","monosaccharide catabolic process",0.3010973056012369,0.5096471993768801,5.7407743136213645,5.087536525907451,-3.7576772647090952,0.9222377691766676,0.60243416),
                     c("GO:0046390","ribose phosphate biosynthetic process",1.4600671821363274,1.9558434254319643,5.405911385852896,5.773199678110848,-3.186251596830003,0.887136849542116,0.68812073),
                     c("GO:0046620","regulation of organ growth",0.02131285513938617,-7.778883935944213,-0.09736398853051287,3.937517892017347,-3.321856896944227,0.7783207985173446,0.63680808),
                     c("GO:0046716","muscle cell cellular homeostasis",0.007270836595651548,4.087347436488592,4.48985423493366,3.4705574852172743,-4.106238237942057,0.9919428007432214,-0),
                     c("GO:0046883","regulation of hormone secretion",0.04974886569780267,-4.380680271261727,-3.7497162015699415,4.305630775996965,-3.954028020939563,0.6943314123874904,0.69827063),
                     c("GO:0048646","anatomical structure formation involved in morphogenesis",0.41648799295422956,2.8217615256752095,-5.9539445506929844,5.228431158637907,-8.200659450546418,0.7983030283082914,0.68222728),
                     c("GO:0048662","negative regulation of smooth muscle cell proliferation",0.005956474123722662,-2.967036932199432,2.73574756865876,3.383994789441733,-3.176853245358756,0.8648005499860926,0.11229583),
                     c("GO:0048769","sarcomerogenesis",0.0014841958250432915,3.926369993962982,-4.269513765449537,2.7810369386211318,-3.2525827537551883,0.7655060791569582,0.69610262),
                     c("GO:0050789","regulation of biological process",26.67823535368239,-5.141615543723029,-1.5380160294482452,7.03498324614998,-3.3267689120150763,0.8113572410450329,0.67600498),
                     c("GO:0050793","regulation of developmental process",1.4293716495756061,-5.734617009378539,-1.2958362398828889,5.763972018350213,-3.455637975765367,0.8519310187670869,0.23698226),
                     c("GO:0050794","regulation of cellular process",24.81324853742716,-5.192820706105049,-1.3531902654658174,7.003509775209796,-3.2078724698868584,0.8045039659008172,0.51240944),
                     c("GO:0050848","regulation of calcium-mediated signaling",0.013291305895910072,-5.260699058824226,3.2873980479987255,3.7324741772811936,-4.7594507517174005,0.867671305829015,0.47134066),
                     c("GO:0050896","response to stimulus",16.94347616187058,5.060394608388128,3.594555715577488,6.83782868559364,-4.247183568811729,1,-0),
                     c("GO:0051128","regulation of cellular component organization",1.4586420587819326,-4.574597927151585,-0.8929144520939848,5.772775571313585,-3.7988706365654488,0.8464485839713456,0.23760414),
                     c("GO:0051147","regulation of muscle cell differentiation",0.022767514729105216,-7.075483738925003,1.2823418184212492,3.966188680956137,-3.3408123528555866,0.8364504130349362,0.64002138),
                     c("GO:0051153","regulation of striated muscle cell differentiation",0.015474525956960488,-7.239019488888379,1.4352508008102443,3.7985125330313516,-3.523365427693123,0.8333358856326215,0.64407125),
                     c("GO:0051239","regulation of multicellular organismal process",0.8600139479948613,-5.4697927429835556,-2.2524900586984864,5.543332844310566,-6.698970004336019,0.8580001599120135,0.20482219),
                     c("GO:0051259","protein complex oligomerization",0.23557116753440024,5.230967917453593,-1.3267122316981554,4.980952778597704,-3.87880139741531,0.95687565380684,0.33684735),
                     c("GO:0051282","regulation of sequestering of calcium ion",0.08941602973826317,-2.8291286642091014,-2.604144668326937,4.560253443542344,-4.346787486224656,0.7524873208805424,0.68961029),
                     c("GO:0051493","regulation of cytoskeleton organization",0.47119894618650526,-4.740134036938125,0.23531173770435918,5.2820326856053805,-5.201349354554731,0.8440381404175918,0.19257431),
                     c("GO:0055074","calcium ion homeostasis",0.27711388251795027,3.9670238341879305,4.610417193947428,5.051488247049932,-3.722835620977585,0.9904112057772171,0.50223633),
                     c("GO:0060047","heart contraction",0.04380715969174213,0.053334362824196034,-8.051354058560928,4.250395603057116,-16.036212172654444,0.8085550186222288,0.61617029),
                     c("GO:0060306","regulation of membrane repolarization",0.01570589313366707,-4.84969926158328,-5.756586908947646,3.8049567998574916,-7.562249437179612,0.8314957002519925,0.62810377),
                     c("GO:0060307","regulation of ventricular cardiac muscle cell membrane repolarization",0.014832112838658167,-4.690131945111699,-5.774798504212278,3.7801011914679115,-6.623423042943488,0.8295076147270678,0.62929165),
                     c("GO:0060973","cell migration involved in heart development",0.0017278697664683093,1.7924133791467824,-6.768640360709556,2.846955325019824,-3.2568564826611626,0.7568131823561873,0.6531505),
                     c("GO:0061061","muscle structure development",0.2045704272086618,2.983351473287219,-6.209158107220858,4.919674183960281,-24.853871964321762,0.8198185117240753,-0),
                     c("GO:0061615","glycolytic process through fructose-6-phosphate",0.05285509311273572,1.040925167961601,5.8005704405441225,4.3319331725032795,-4.616184634019569,0.8612936484429584,0.52189436),
                     c("GO:0065007","biological regulation",27.969320810822822,5.017293182263625,2.6828891650502897,7.055508065118343,-3.1573090258726473,1,-0),
                     c("GO:0065008","regulation of biological quality",3.028571729559731,-4.744221563967171,-1.9467527987854394,6.0900643235623235,-4.987162775294828,0.8488083800066717,0.22399414),
                     c("GO:0070252","actin-mediated cell contraction",0.020682748785802285,4.880872920708508,-3.1451131290022847,3.924486043733915,-12.349692476868064,0.9192886419424964,0.63052256),
                     c("GO:0070885","negative regulation of calcineurin-NFAT signaling cascade",0.0026533384732946408,-4.593282405249216,3.6954162550235354,3.0330214446829107,-6.0762380391713,0.8465784584645704,0.34417553),
                     c("GO:0071415","cellular response to purine-containing compound",0.0009451595303758273,-1.9266936623729818,6.261645738897356,2.5854607295085006,-5.34008379993015,0.9742814322740034,0.18899936),
                     c("GO:0072359","circulatory system development",0.32041400350329285,1.8082214141531001,-6.978099363844881,5.114540931085493,-8.835647144215564,0.7257315017703208,0.6271619),
                     c("GO:0086001","cardiac muscle cell action potential",0.014347226308751817,-4.722337996051299,-6.135649441761177,3.765668554759014,-7.863279432843593,0.8309473496291302,0.11208298),
                     c("GO:0086036","regulation of cardiac muscle cell membrane potential",0.00166879729581982,-4.5083283935209515,-6.813651613344267,2.8318697742805017,-4.0746879085003505,0.8542239113538314,0.55716408),
                     c("GO:0090066","regulation of anatomical structure size",0.38626011745280886,-4.963347406314929,-4.500413935168436,5.195708742179063,-3.936813057124898,0.8192040732114312,0.4708064),
                     c("GO:0097435","supramolecular fiber organization",0.789545413217101,5.8697877627685076,-1.3584253598820395,5.5062045753034194,-12.703334809738468,0.9626653370011518,0.3181356),
                     c("GO:0097623","potassium ion export across plasma membrane",0.008747648361863777,3.4706662764999177,1.5100846833381123,3.550839605065785,-3.2594170606181403,0.989661209407211,0.59171161),
                     c("GO:0140115","export across plasma membrane",0.15693586368948637,3.3059989303661,1.3239141861504453,4.80455511972865,-4.071092309756048,0.9884306975942023,0.01325593),
                     c("GO:1901880","negative regulation of protein depolymerization",0.09490730815562899,-4.298752803748315,1.4201091125344731,4.586137025230793,-3.093022612968563,0.8337137410437306,0.67642628),
                     c("GO:1903115","regulation of actin filament-based movement",0.004602730004694785,-0.15272284468804792,-0.6231882556146767,3.27207378750001,-3.224133346511831,0.8918577245357027,0.11064255),
                     c("GO:1903814","regulation of collecting lymphatic vessel constriction",4.922705887374102E-06,-6.119025153467328,-6.161497465930297,0.47712125471966244,-3.25276691892186,0.8229651691745878,0.53287473),
                     c("GO:1904062","regulation of monoatomic cation transmembrane transport",0.1504083556828283,-3.1789701669383925,-2.522311029222786,4.786105176895488,-8.903089986991944,0.7056323292483629,0.13814733),
                     c("GO:1904753","negative regulation of vascular associated smooth muscle cell migration",0.0005464203534985252,-2.1986092652536096,3.449632874455222,2.3483048630481607,-4.3260580013659125,0.8718194286627357,0.62053967),
                     c("GO:2000145","regulation of cell motility",0.3916726325759767,-3.614546717764385,-0.11743203009260395,5.201752062778252,-4.716698771296451,0.8443560706302523,0.18910017))


revigo.data_2 <- rbind(c("GO:0001725","stress fiber",0.06008236809071794,6.803249517776746,2.159492763552769,4.381638446726157,-9.176525770829699,0.6664047305574621,0.66519672),
                     c("GO:0005856","cytoskeleton",3.1141546278807093,7.44461247563631,1.1012630002452035,6.096213889361101,-8.943095148663527,0.7140177438011025,0.49045169),
                     c("GO:0005859","muscle myosin complex",0.0006188399072388922,5.328500424048609,1.4087552921651532,2.3961993470957363,-6.5934598195660445,0.5375432331303068,0.67653506),
                     c("GO:0005862","muscle thin filament tropomyosin",2.4953222066084366E-05,5.4824279998682615,0.6670461297602915,1.0413926851582251,-3.3751694509421903,0.5870801202588505,0.58090183),
                     c("GO:0005863","striated muscle myosin thick filament",0.0002520275428674521,6.36713147932564,-1.1568086912025284,2.0086001717619175,-4.876148359032914,0.6262802946855078,0.65054002),
                     c("GO:0005886","plasma membrane",16.678688713171073,-3.15490123482533,5.033539515122577,6.8250353366409255,-4.341988603342887,0.9429604568494633,0.31779959),
                     c("GO:0014704","intercalated disc",0.014193392711188788,-0.3269780274679795,-5.8024639196837375,3.7550359337677714,-9.707743928643524,0.7951884383132329,0.60552702),
                     c("GO:0015629","actin cytoskeleton",0.7637881835763632,6.919746044261761,1.564442544395725,5.485847722524024,-16.13608262304214,0.6753014568138466,0.36298147),
                     c("GO:0016010","dystrophin-associated glycoprotein complex",0.028451663799749394,-0.6714249575026296,5.894419390114396,4.057019124322766,-3.5977712983596875,0.8701135726294176,0.46246365),
                     c("GO:0016460","myosin II complex",0.07494201183107119,5.130778529245161,3.0629875808241063,4.477613176429474,-6.812479279163537,0.6221932249278012,0.67546024),
                     c("GO:0016528","sarcoplasm",0.05873239877694277,-5.592926439377034,-0.31187214173753247,4.371769558513179,-3.74865140575919,0.9513975161173657,0.07342577),
                     c("GO:0030017","sarcomere",0.16626082330411354,6.241144672092539,-0.2972598411362188,4.823669813268136,-33.20273245916928,0.5001604342516223,0),
                     c("GO:0030054","cell junction",1.8830924078948503,4.404703284453421,-6.835340417602166,5.87774557638143,-5.431798275933005,0.9999356900717639,4.293E-05),
                     c("GO:0030314","junctional membrane complex",0.007817844473304233,3.1334761806935854,7.718020132212226,3.496098992132571,-4.357535479757878,0.8798177656653938,0.18553217),
                     c("GO:0030315","T-tubule",0.01728509692517664,-3.8496187830368203,4.34545456946291,3.84060787900929,-4.482804102050026,0.9632540181144719,0.19190016),
                     c("GO:0030425","dendrite",0.3835983968552951,-3.7907930809254977,-4.326811662943234,5.186752977050775,-3.332137229181728,0.9673519306671048,3.581E-05),
                     c("GO:0032279","asymmetric synapse",0.18297948208839007,0.00366597520296921,-6.096671266959751,4.865281684995611,-4.573488738635425,0.7546065216151935,0.53367494),
                     c("GO:0032432","actin filament bundle",0.08030695457527932,6.4112043022084375,1.9238874266623738,4.507640019563205,-8.61261017366127,0.6695928091284603,0.56011642),
                     c("GO:0032982","myosin filament",0.04267749569962409,5.91617173099078,-3.15400551970002,4.233097687864703,-4.2335871528876,0.6856864914259838,0.60980044),
                     c("GO:0034703","cation channel complex",0.6343358581419307,1.1080523295607987,6.312244971294455,5.405194339042305,-5.347753658996677,0.8465898997359792,0.0583811),
                     c("GO:0042383","sarcolemma",0.10179417409638457,-3.3944787316104463,5.629607406224533,4.610606937465461,-12.028724151261894,0.9586363470614938,2.583E-05),
                     c("GO:0042599","lamellar body",0.0027124152385833706,-2.3722288690822695,-1.966774944291344,3.036628895362161,-3.762707662432541,0.9227195353324413,0.08004435),
                     c("GO:0042641","actomyosin",0.0625552323974669,6.491581140227685,2.218828200712562,4.399154333958217,-8.27490547891853,0.6743454327165566,0.66704654),
                     c("GO:0042995","cell projection",2.454214268576769,-5.6032390733422455,1.9128029795107933,5.9927862867668535,-3.4460793684788302,0.999933941737404,4.469E-05),
                     c("GO:0043292","contractile muscle fiber",0.1962146710722412,6.615131982736691,-0.2791638235717606,4.895610368172215,-34.6252516539899,0.5713749130857955,0.68343832),
                     c("GO:0044291","cell-cell contact zone",0.02478853080044821,0.4248025282557036,-6.1105442298042805,3.997167871445834,-10.347753658996677,0.7882446087950364,2.335E-05),
                     c("GO:0071944","cell periphery",18.431522806558757,-3.0976109172220156,0.560078613413288,6.868434645240301,-5.350665141287858,0.9999168396270548,6.358E-05),
                     c("GO:0090665","glycoprotein complex",0.0287785510088151,1.3257506203254086,7.574874913734424,4.061979947074878,-3.5977712983596875,0.9333630292220478,0.20383064),
                     c("GO:0097512","cardiac myofibril",0.0005988773295860248,5.984573380636086,-0.983110091114557,2.3820170425748683,-11.453457336521868,0.6223686729350679,0.68109929),
                     c("GO:0098984","neuron to neuron synapse",0.18595640148087392,0.11672865685635982,-5.80724362525155,4.872290329547083,-4.543633966870957,0.7633822502203861,0.69550658),
                     c("GO:0099080","supramolecular complex",1.7908628038163958,-4.985497209696189,-2.451035187079568,5.855936289863483,-26.429457060118104,0.9999360109819824,3.293E-05),
                     c("GO:1990351","transporter complex",2.0712297209843,2.132033936507174,6.618511999489699,5.919102161041849,-3.879896809513499,0.909911296279212,0.30135753))

revigo.data_3 <- rbind(c("GO:0003779","actin binding",0.7436112368550806,-6.234706682965,0.58084494878162,5.532345363784474,-12.821023052706831,0.5945608347104586,0.50318841),
                     c("GO:0005004","GPI-linked ephrin receptor activity",0.00040817223680983226,2.279428565516859,-1.5491746485657674,2.27415784926368,-3.2461099031293466,0.9799961640258203,0.00818148),
                     c("GO:0005179","hormone activity",0.12371765979883045,-4.995691162654582,-1.0707994182736693,4.753437503725226,-3.6246400411561996,0.6454604368591427,0.57060577),
                     c("GO:0005198","structural molecule activity",3.0956240814896883,4.029729606819689,-4.031588629342251,6.151746667822227,-7.667561540084395,1,-0),
                     c("GO:0005200","structural constituent of cytoskeleton",0.14988652047857953,0.05485862474763788,6.959445868830417,4.836767047394205,-4.379863945026242,0.9409196177136391,0.44086029),
                     c("GO:0005488","binding",58.52264612691586,5.1383183929395235,-0.881030490920006,7.428322116042559,-3.4799699947543847,1,-0),
                     c("GO:0005515","protein binding",8.27054301946569,4.63992372054119,1.774768946575108,6.578532284734752,-10.838631997765026,0.9538052696150507,0.06654784),
                     c("GO:0005539","glycosaminoglycan binding",0.16800107338384584,3.7651939891302195,3.963318679800449,4.886315844136417,-3.1193007107812982,0.9661255137886532,0.0420273),
                     c("GO:0008092","cytoskeletal protein binding",1.5207580793055957,-6.5193245025054205,-0.2983837523645333,5.843058900884321,-18.02826040911222,0.6790740866308421,0),
                     c("GO:0008307","structural constituent of muscle",0.015464707474853807,-0.5246910812149652,7.242635874073894,3.8504011479971583,-13.962573502059376,0.9513464730438114,-0),
                     c("GO:0008324","monoatomic cation transmembrane transporter activity",2.870764833792156,-1.1022296565915843,-7.359501550413587,6.118996092927394,-3.1552427045817253,0.7034046149837607,0.6956531),
                     c("GO:0022843","voltage-gated monoatomic cation channel activity",0.28304889294412733,-0.7933142932281428,-7.260968544034136,5.112862954813035,-4.801342913045577,0.6530921979781368,-0),
                     c("GO:0030506","ankyrin binding",0.003289388026055708,-6.751464053658378,1.9145268586613389,3.178401341533755,-3.3746793160456168,0.6981614931245017,0.54156632),
                     c("GO:0031432","titin binding",0.0017571050835931284,-5.923369710520075,2.2347639898746143,2.906335041805091,-3.5964419716301745,0.7070665948218637,0.52014258),
                     c("GO:0032036","myosin heavy chain binding",0.01739424895795483,-7.044236517677472,0.9478997344048669,3.9014583213961123,-3.7418442034344834,0.6718116672540403,0.60809031),
                     c("GO:0035256","G protein-coupled glutamate receptor binding",0.005675122009120663,-6.122235130291486,-1.9876965115032137,3.4151403521958725,-3.2522107763363772,0.7049268091111228,0.67089612),
                     c("GO:0042805","actinin binding",0.02319815258189785,-6.443387586161278,1.2025336987461344,4.02649240705284,-10.536107011014092,0.6668154128055556,0.62128442),
                     c("GO:0044325","transmembrane transporter binding",0.06262278863141224,-4.729970453268283,0.09970261795154277,4.457745685468492,-7.649751981665837,0.7380450579338944,0.39549856),
                     c("GO:0051371","muscle alpha-actinin binding",0.01762780205602249,-5.736373817098115,1.371670321478236,3.9072500828813284,-11.43297363384094,0.6457925455447343,0.6086891),
                     c("GO:0051428","peptide hormone receptor binding",0.01288034422146963,-6.741265358490199,-1.5894164759937917,3.7709992051639407,-4.505845405981558,0.7022202432547745,0.347911),
                     c("GO:0071855","neuropeptide receptor binding",0.029322918873279608,-6.083572077834979,-1.4457556752855978,4.128237670769187,-3.6037117101307143,0.6799389111452252,0.60226019),
                     c("GO:0097493","structural molecule activity conferring elasticity",0.0005565985047406804,0.7463549993043072,6.546226493558778,2.4082399653118496,-4.048662481204082,0.9563109078601626,0.32267231),
                     c("GO:0099580","monoatomic ion antiporter activity involved in regulation of postsynaptic membrane potential",2.182739234277178E-06,0.9699669325229828,-6.410795089494994,0.3010299956639812,-3.2497925528209737,0.9074531528728198,0.17759385))



# Combine both data frames and label them with categories
revigo.data_1 <- data.frame(revigo.data_1)
revigo.data_2 <- data.frame(revigo.data_2)
revigo.data_3 <- data.frame(revigo.data_3)

names(revigo.data_1) <- revigo.names
names(revigo.data_2) <- revigo.names
names(revigo.data_3) <- revigo.names
revigo.data_1$Category <- "Biological process"
revigo.data_2$Category <- "Cellular component"
revigo.data_3$Category <- "Molecular function"

# Combine both data frames
combined_data <- rbind(revigo.data_1, revigo.data_2, revigo.data_3)

# Clean the data (convert to numeric)
combined_data$plot_X <- as.numeric(as.character(combined_data$plot_X))
combined_data$plot_Y <- as.numeric(as.character(combined_data$plot_Y))
combined_data$log_size <- as.numeric(as.character(combined_data$log_size))
combined_data$value <- as.numeric(as.character(combined_data$value))
combined_data$frequency <- as.numeric(as.character(combined_data$frequency))
combined_data$uniqueness <- as.numeric(as.character(combined_data$uniqueness))
combined_data$dispensability <- as.numeric(as.character(combined_data$dispensability))

# Make three faceted plots. Note we are only displaying the "cluster represantative" term 

# Filter data for specific descriptions (cluster representative)
filtered_data <- subset(combined_data, description %in% c("GPI-linked ephrin receptor activity", "hormone activity", "structural molecule activity", "binding", "glycosaminoglycan binding", "cytoskeletal protein binding", "structural constituent of muscle", "voltage-gated monoatomic cation channel activity", "sarcomere", "cell junction", "dendrite", "sarcolemma", "cell projection", "cell-cell contact zone", "cell periphery", "supramolecular complex", "muscle system process", "pyruvate metabolic process", "cell communication", "regulation of heart contraction", "response to wounding", "actin filament-based process", "actomyosin structure organization", "regulation of transporter activity", "multicellular organismal process", "developmental process", "locomotion", "homeostatic process", "muscle cell cellular homeostasis", "regulation of hormone secretion", "response to stimulus", "muscle structure development", "biological regulation", "export across plasma membrane"))

# Ensure 'Category' is a factor with levels ordered as desired
combined_data$Category <- factor(combined_data$Category, 
                                 levels = c("Biological process", "Cellular component", "Molecular function"))
p_faceted <- ggplot(data = combined_data) +
  geom_point(aes(plot_X, plot_Y, colour = Category, size = value), alpha = 0.6) +
  scale_colour_manual(values = c("Cellular component" = "#fbb4ae", 
                                 "Biological process" = "#b3cde3", 
                                 "Molecular function" = "#ccebc5")) +
  # Reduce number of size legend breaks
  scale_size_continuous(range = c(3, 15), breaks = c(-34.6, -27.8, -21, -14.2, -7.4)) +  
  geom_point(aes(plot_X, plot_Y, size = value), shape = 21, 
             fill = "transparent", colour = I(adjustcolor("black", alpha.f = 0.6))) +
  # Add labels only for filtered data with ggrepel
  geom_label_repel(data = filtered_data, aes(plot_X, plot_Y, label = description), 
                   fill = "white", colour = I(adjustcolor("black", alpha.f = 0.85)), 
                   size = 5, nudge_y = -1, 
                   max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3) +
  labs(y = "semantic space x", x = "semantic space y", size = "log(adjusted p value)") +
  theme_bw() +
  theme(legend.key = element_blank(),
        legend.text = element_text(size = 18),  # Increase the size of legend text
        legend.title = element_text(size = 20), # Increase the size of legend titles
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        strip.text = element_text(size = 22),
        plot.margin = margin(10, 10, 10, 10),         # Add margins around the plot
        panel.spacing = unit(1.5, "lines"),
        axis.text.x = element_text(margin = margin(t = 5))  # Add extra space at the top of x-axis labels
  ) +
  xlim(min(combined_data$plot_X) - 1, max(combined_data$plot_X) + 1) +  # Add buffer
  ylim(min(combined_data$plot_Y) - 1, max(combined_data$plot_Y) + 1) +  # Add buffer
  facet_wrap(~Category, scales = "free", nrow = 1) +  # Two rows for facets
  guides(colour = guide_legend(override.aes = list(size = 16, shape = 21, fill = c("#fbb4ae", "#b3cde3", "#ccebc5"), colour = "black")))  # Add black outline to legend

p_faceted 

ggsave("semantic_full_freeze_up_over.pdf", p_faceted, width = 32, height = 8, units = "in")



