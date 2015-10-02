! This file is F-compatible, except for upper/lower case conventions.
!--------------------------------------------------------------------
Module QuadratureRule

USE Precision_Model, ONLY: stnd

Implicit NONE

PRIVATE

PUBLIC :: Dqk_drv
PRIVATE :: Dqknn

CONTAINS
SUBROUTINE Dqk_drv(KEY,VER,NUMFUN,Integrand,BASVAL,RGNERR,NUM)
!***
! Driver routine voor Dqknn. This routine selects the appropriate data
! for subroutine Dqknn, according to the value of key, i.e. the 
! desired integration rule and error estimatation routine.
!
!*** 
! For the meaning of the parameters, see Dqknn
!
      INTEGER, INTENT(IN) :: NUMFUN,KEY
      INTEGER, INTENT(OUT) :: NUM
      REAL(kind=stnd), DIMENSION(:,:),INTENT(IN) :: VER
      REAL(kind=stnd), DIMENSION(:),INTENT(OUT) :: BASVAL, RGNERR
      INTERFACE 
         FUNCTION Integrand(NUMFUN,X) RESULT(Value)
            USE Precision_Model
            INTEGER, INTENT(IN) :: NUMFUN
            REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
            REAL(kind=stnd), DIMENSION(NUMFUN) :: Value
         END FUNCTION Integrand
      END INTERFACE
!
      REAL(kind=stnd), DIMENSION(1:15) :: wg
      REAL(kind=stnd), DIMENSION(1:31) :: xgk, wgk
!
!***first executable statement  dqk_drv
!
   SELECT CASE (key)

   CASE(:1)
    wg(1:4) = (/ 0.129484966168869693270611432679082_stnd,     &
                 0.279705391489276667901467771423780_stnd,     &
                 0.381830050505118944950369775488975_stnd,     &
                 0.417959183673469387755102040816327_stnd     /)
    xgk(1:8) = (/ 0.991455371120812639206854697526329_stnd,    &
                  0.949107912342758524526189684047851_stnd,    &
                  0.864864423359769072789712788640926_stnd,    &
                  0.741531185599394439863864773280788_stnd,    &
                  0.586087235467691130294144838258730_stnd,    &
                  0.405845151377397166906606412076961_stnd,    &
                  0.207784955007898467600689403773245_stnd,    &
                  0.000000000000000000000000000000000_stnd    /)
    wgk(1:8) = (/ 0.022935322010529224963732008058970_stnd,    &
                  0.063092092629978553290700663189204_stnd,    &
                  0.104790010322250183839876322541518_stnd,    &
                  0.140653259715525918745189590510238_stnd,    &
                  0.169004726639267902826583426598550_stnd,    &
                  0.190350578064785409913256402421014_stnd,    &
                  0.204432940075298892414161999234649_stnd,    &
                  0.209482141084727828012999174891714_stnd    /)
    NUM = 15
    call Dqknn(KEY,VER,NUMFUN,Integrand,BASVAL,RGNERR,  &
               wg(1:4),xgk(1:8),wgk(1:8))

   CASE (2)
    wg(1:5) = (/ 0.066671344308688137593568809893332_stnd,     &
                 0.149451349150580593145776339657697_stnd,     &
                 0.219086362515982043995534934228163_stnd,     &
                 0.269266719309996355091226921569469_stnd,     &
                 0.295524224714752870173892994651338_stnd     /)
     xgk(1:11) = (/ 0.995657163025808080735527280689003_stnd,  &
                    0.973906528517171720077964012084452_stnd,  &
                    0.930157491355708226001207180059508_stnd,  &
                    0.865063366688984510732096688423493_stnd,  &
                    0.780817726586416897063717578345042_stnd,  &
                    0.679409568299024406234327365114874_stnd,  &
                    0.562757134668604683339000099272694_stnd,  &
                    0.433395394129247190799265943165784_stnd,  &
                    0.294392862701460198131126603103866_stnd,  &
                    0.148874338981631210884826001129720_stnd,  &
                    0.000000000000000000000000000000000_stnd  /)
     wgk(1:11) = (/ 0.011694638867371874278064396062192_stnd,  &
                    0.032558162307964727478818972459390_stnd,  &
                    0.054755896574351996031381300244580_stnd,  &
                    0.075039674810919952767043140916190_stnd,  &
                    0.093125454583697605535065465083366_stnd,  &
                    0.109387158802297641899210590325805_stnd,  &
                    0.123491976262065851077958109831074_stnd,  &
                    0.134709217311473325928054001771707_stnd,  &
                    0.142775938577060080797094273138717_stnd,  &
                    0.147739104901338491374841515972068_stnd,  &
                    0.149445554002916905664936468389821_stnd  /)
    NUM = 21
    call Dqknn(KEY,VER,NUMFUN,Integrand,BASVAL,RGNERR,  &
               wg(1:5),xgk(1:11),wgk(1:11))

   CASE (3)
    wg(1:8) = (/ 0.030753241996117268354628393577204_stnd,  &
                 0.070366047488108124709267416450667_stnd,  &
                 0.107159220467171935011869546685869_stnd,  &
                 0.139570677926154314447804794511028_stnd,  &
                 0.166269205816993933553200860481209_stnd,  &
                 0.186161000015562211026800561866423_stnd,  &
                 0.198431485327111576456118326443839_stnd,  &
                 0.202578241925561272880620199967519_stnd  /)
    xgk(1:16) = (/ 0.998002298693397060285172840152271_stnd,  &
                   0.987992518020485428489565718586613_stnd,  &
                   0.967739075679139134257347978784337_stnd,  &
                   0.937273392400705904307758947710209_stnd,  &
                   0.897264532344081900882509656454496_stnd,  &
                   0.848206583410427216200648320774217_stnd,  &
                   0.790418501442465932967649294817947_stnd,  &
                   0.724417731360170047416186054613938_stnd,  &
                   0.650996741297416970533735895313275_stnd,  &
                   0.570972172608538847537226737253911_stnd,  &
                   0.485081863640239680693655740232351_stnd,  &
                   0.394151347077563369897207370981045_stnd,  &
                   0.299180007153168812166780024266389_stnd,  &
                   0.201194093997434522300628303394596_stnd,  &
                   0.101142066918717499027074231447392_stnd,  &
                   0.000000000000000000000000000000000_stnd  /)
    wgk(1:16) = (/ 0.005377479872923348987792051430128_stnd,  &
                   0.015007947329316122538374763075807_stnd,  &
                   0.025460847326715320186874001019653_stnd,  &
                   0.035346360791375846222037948478360_stnd,  &
                   0.044589751324764876608227299373280_stnd,  &
                   0.053481524690928087265343147239430_stnd,  &
                   0.062009567800670640285139230960803_stnd,  &
                   0.069854121318728258709520077099147_stnd,  &
                   0.076849680757720378894432777482659_stnd,  &
                   0.083080502823133021038289247286104_stnd,  &
                   0.088564443056211770647275443693774_stnd,  &
                   0.093126598170825321225486872747346_stnd,  &
                   0.096642726983623678505179907627589_stnd,  &
                   0.099173598721791959332393173484603_stnd,  &
                   0.100769845523875595044946662617570_stnd,  &
                   0.101330007014791549017374792767493_stnd  /)
    NUM = 31
    call Dqknn(KEY,VER,NUMFUN,Integrand,BASVAL,RGNERR,  &
               wg(1:8),xgk(1:16),wgk(1:16))

   CASE (4)
    wg(1:10) = (/ 0.017614007139152118311861962351853_stnd,  &
                  0.040601429800386941331039952274932_stnd,  &
                  0.062672048334109063569506535187042_stnd,  &
                  0.083276741576704748724758143222046_stnd,  &
                  0.101930119817240435036750135480350_stnd,  &
                  0.118194531961518417312377377711382_stnd,  &
                  0.131688638449176626898494499748163_stnd,  &
                  0.142096109318382051329298325067165_stnd,  &
                  0.149172986472603746787828737001969_stnd,  &
                  0.152753387130725850698084331955098_stnd  /)
    xgk(1:21) = (/ 0.998859031588277663838315576545863_stnd,  &
                   0.993128599185094924786122388471320_stnd,  &
                   0.981507877450250259193342994720217_stnd,  &
                   0.963971927277913791267666131197277_stnd,  &
                   0.940822633831754753519982722212443_stnd,  &
                   0.912234428251325905867752441203298_stnd,  &
                   0.878276811252281976077442995113078_stnd,  &
                   0.839116971822218823394529061701521_stnd,  &
                   0.795041428837551198350638833272788_stnd,  &
                   0.746331906460150792614305070355642_stnd,  &
                   0.693237656334751384805490711845932_stnd,  &
                   0.636053680726515025452836696226286_stnd,  &
                   0.575140446819710315342946036586425_stnd,  &
                   0.510867001950827098004364050955251_stnd,  &
                   0.443593175238725103199992213492640_stnd,  &
                   0.373706088715419560672548177024927_stnd,  &
                   0.301627868114913004320555356858592_stnd,  &
                   0.227785851141645078080496195368575_stnd,  &
                   0.152605465240922675505220241022678_stnd,  &
                   0.076526521133497333754640409398838_stnd,  &
                   0.000000000000000000000000000000000_stnd  /)
    wgk(1:21) = (/ 0.003073583718520531501218293246031_stnd,  &
                   0.008600269855642942198661787950102_stnd,  &
                   0.014626169256971252983787960308868_stnd,  &
                   0.020388373461266523598010231432755_stnd,  &
                   0.025882133604951158834505067096153_stnd,  &
                   0.031287306777032798958543119323801_stnd,  &
                   0.036600169758200798030557240707211_stnd,  &
                   0.041668873327973686263788305936895_stnd,  &
                   0.046434821867497674720231880926108_stnd,  &
                   0.050944573923728691932707670050345_stnd,  &
                   0.055195105348285994744832372419777_stnd,  &
                   0.059111400880639572374967220648594_stnd,  &
                   0.062653237554781168025870122174255_stnd,  &
                   0.065834597133618422111563556969398_stnd,  &
                   0.068648672928521619345623411885368_stnd,  &
                   0.071054423553444068305790361723210_stnd,  &
                   0.073030690332786667495189417658913_stnd,  &
                   0.074582875400499188986581418362488_stnd,  &
                   0.075704497684556674659542775376617_stnd,  &
                   0.076377867672080736705502835038061_stnd,  &
                   0.076600711917999656445049901530102_stnd  /)
    NUM = 41
    call Dqknn(KEY,VER,NUMFUN,Integrand,BASVAL,RGNERR,  &
               wg(1:10),xgk(1:21),wgk(1:21))

   CASE (5)
    wg(1:13) = (/ 0.011393798501026287947902964113235_stnd,  &
                  0.026354986615032137261901815295299_stnd,  &
                  0.040939156701306312655623487711646_stnd,  &
                  0.054904695975835191925936891540473_stnd,  &
                  0.068038333812356917207187185656708_stnd,  &
                  0.080140700335001018013234959669111_stnd,  &
                  0.091028261982963649811497220702892_stnd,  &
                  0.100535949067050644202206890392686_stnd,  &
                  0.108519624474263653116093957050117_stnd,  &
                  0.114858259145711648339325545869556_stnd,  &
                  0.119455763535784772228178126512901_stnd,  &
                  0.122242442990310041688959518945852_stnd,  &
                  0.123176053726715451203902873079050_stnd  /)
    xgk(1:26) = (/ 0.999262104992609834193457486540341_stnd,  &
                   0.995556969790498097908784946893902_stnd,  &
                   0.988035794534077247637331014577406_stnd,  &
                   0.976663921459517511498315386479594_stnd,  &
                   0.961614986425842512418130033660167_stnd,  &
                   0.942974571228974339414011169658471_stnd,  &
                   0.920747115281701561746346084546331_stnd,  &
                   0.894991997878275368851042006782805_stnd,  &
                   0.865847065293275595448996969588340_stnd,  &
                   0.833442628760834001421021108693570_stnd,  &
                   0.797873797998500059410410904994307_stnd,  &
                   0.759259263037357630577282865204361_stnd,  &
                   0.717766406813084388186654079773298_stnd,  &
                   0.673566368473468364485120633247622_stnd,  &
                   0.626810099010317412788122681624518_stnd,  &
                   0.577662930241222967723689841612654_stnd,  &
                   0.526325284334719182599623778158010_stnd,  &
                   0.473002731445714960522182115009192_stnd,  &
                   0.417885382193037748851814394594572_stnd,  &
                   0.361172305809387837735821730127641_stnd,  &
                   0.303089538931107830167478909980339_stnd,  &
                   0.243866883720988432045190362797452_stnd,  &
                   0.183718939421048892015969888759528_stnd,  &
                   0.122864692610710396387359818808037_stnd,  &
                   0.061544483005685078886546392366797_stnd,  &
                   0.000000000000000000000000000000000_stnd  /)
    wgk(1:26) = (/ 0.001987383892330315926507851882843_stnd,  &
                   0.005561932135356713758040236901066_stnd,  &
                   0.009473973386174151607207710523655_stnd,  &
                   0.013236229195571674813656405846976_stnd,  &
                   0.016847817709128298231516667536336_stnd,  &
                   0.020435371145882835456568292235939_stnd,  &
                   0.024009945606953216220092489164881_stnd,  &
                   0.027475317587851737802948455517811_stnd,  &
                   0.030792300167387488891109020215229_stnd,  &
                   0.034002130274329337836748795229551_stnd,  &
                   0.037116271483415543560330625367620_stnd,  &
                   0.040083825504032382074839284467076_stnd,  &
                   0.042872845020170049476895792439495_stnd,  &
                   0.045502913049921788909870584752660_stnd,  &
                   0.047982537138836713906392255756915_stnd,  &
                   0.050277679080715671963325259433440_stnd,  &
                   0.052362885806407475864366712137873_stnd,  &
                   0.054251129888545490144543370459876_stnd,  &
                   0.055950811220412317308240686382747_stnd,  &
                   0.057437116361567832853582693939506_stnd,  &
                   0.058689680022394207961974175856788_stnd,  &
                   0.059720340324174059979099291932562_stnd,  &
                   0.060539455376045862945360267517565_stnd,  &
                   0.061128509717053048305859030416293_stnd,  &
                   0.061471189871425316661544131965264_stnd,  &
!       note: wgk (26) was calculated from the values of wgk(1..25)
                   0.061580818067832935078759824240066_stnd  /)
    NUM = 51
    call Dqknn(KEY,VER,NUMFUN,Integrand,BASVAL,RGNERR,  &
               wg(1:13),xgk(1:26),wgk(1:26))

   CASE (6:)
    wg(1:15) = (/ 0.007968192496166605615465883474674_stnd,  &
                  0.018466468311090959142302131912047_stnd,  &
                  0.028784707883323369349719179611292_stnd,  &
                  0.038799192569627049596801936446348_stnd,  &
                  0.048402672830594052902938140422808_stnd,  &
                  0.057493156217619066481721689402056_stnd,  &
                  0.065974229882180495128128515115962_stnd,  &
                  0.073755974737705206268243850022191_stnd,  &
                  0.080755895229420215354694938460530_stnd,  &
                  0.086899787201082979802387530715126_stnd,  &
                  0.092122522237786128717632707087619_stnd,  &
                  0.096368737174644259639468626351810_stnd,  &
                  0.099593420586795267062780282103569_stnd,  &
                  0.101762389748405504596428952168554_stnd,  &
                  0.102852652893558840341285636705415_stnd  /)
    xgk(1:31) = (/ 0.999484410050490637571325895705811_stnd,  &
                   0.996893484074649540271630050918695_stnd,  &
                   0.991630996870404594858628366109486_stnd,  &
                   0.983668123279747209970032581605663_stnd,  &
                   0.973116322501126268374693868423707_stnd,  &
                   0.960021864968307512216871025581798_stnd,  &
                   0.944374444748559979415831324037439_stnd,  &
                   0.926200047429274325879324277080474_stnd,  &
                   0.905573307699907798546522558925958_stnd,  &
                   0.882560535792052681543116462530226_stnd,  &
                   0.857205233546061098958658510658944_stnd,  &
                   0.829565762382768397442898119732502_stnd,  &
                   0.799727835821839083013668942322683_stnd,  &
                   0.767777432104826194917977340974503_stnd,  &
                   0.733790062453226804726171131369528_stnd,  &
                   0.697850494793315796932292388026640_stnd,  &
                   0.660061064126626961370053668149271_stnd,  &
                   0.620526182989242861140477556431189_stnd,  &
                   0.579345235826361691756024932172540_stnd,  &
                   0.536624148142019899264169793311073_stnd,  &
                   0.492480467861778574993693061207709_stnd,  &
                   0.447033769538089176780609900322854_stnd,  &
                   0.400401254830394392535476211542661_stnd,  &
                   0.352704725530878113471037207089374_stnd,  &
                   0.304073202273625077372677107199257_stnd,  &
                   0.254636926167889846439805129817805_stnd,  &
                   0.204525116682309891438957671002025_stnd,  &
                   0.153869913608583546963794672743256_stnd,  &
                   0.102806937966737030147096751318001_stnd,  &
                   0.051471842555317695833025213166723_stnd,  &
                   0.000000000000000000000000000000000_stnd  /)
    wgk(1:31) = (/ 0.001389013698677007624551591226760_stnd,  &
                   0.003890461127099884051267201844516_stnd,  &
                   0.006630703915931292173319826369750_stnd,  &
                   0.009273279659517763428441146892024_stnd,  &
                   0.011823015253496341742232898853251_stnd,  &
                   0.014369729507045804812451432443580_stnd,  &
                   0.016920889189053272627572289420322_stnd,  &
                   0.019414141193942381173408951050128_stnd,  &
                   0.021828035821609192297167485738339_stnd,  &
                   0.024191162078080601365686370725232_stnd,  &
                   0.026509954882333101610601709335075_stnd,  &
                   0.028754048765041292843978785354334_stnd,  &
                   0.030907257562387762472884252943092_stnd,  &
                   0.032981447057483726031814191016854_stnd,  &
                   0.034979338028060024137499670731468_stnd,  &
                   0.036882364651821229223911065617136_stnd,  &
                   0.038678945624727592950348651532281_stnd,  &
                   0.040374538951535959111995279752468_stnd,  &
                   0.041969810215164246147147541285970_stnd,  &
                   0.043452539701356069316831728117073_stnd,  &
                   0.044814800133162663192355551616723_stnd,  &
                   0.046059238271006988116271735559374_stnd,  &
                   0.047185546569299153945261478181099_stnd,  &
                   0.048185861757087129140779492298305_stnd,  &
                   0.049055434555029778887528165367238_stnd,  &
                   0.049795683427074206357811569379942_stnd,  &
                   0.050405921402782346840893085653585_stnd,  &
                   0.050881795898749606492297473049805_stnd,  &
                   0.051221547849258772170656282604944_stnd,  &
                   0.051426128537459025933862879215781_stnd,  &
                   0.051494729429451567558340433647099_stnd  /)
    NUM = 61
    call Dqknn(KEY,VER,NUMFUN,Integrand,BASVAL,RGNERR,  &
               wg(1:15),xgk(1:31),wgk(1:31))

   END SELECT

   RETURN
END SUBROUTINE Dqk_drv
!------------------------------------------------------------------------
SUBROUTINE Dqknn(KEY,VER,NUMFUN,Integrand,BASVAL,RGNERR,  &
                 wg,xgk,wgk)
!***begin prologue  dqk15
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***revision date  980326   (yymmdd)  (Integration in Cubpack)
!***revision date  990525   (yymmdd)  (F conversion)
!***category no.  h2a1a2
!***keywords  15-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div - k.u.leuven
!***purpose  to compute i = integral of f over (a,b), with error
!                           estimate
!***description
!
!           integration rules
!           standard fortran subroutine
!           double precision version
!
!           parameters
!   ON ENTRY
!
!   Integrand Externally declared subroutine for computing
!             all components of the integrand at the given
!             evaluation point.
!             It must have parameters (DIMENS,X,NUMFUN,FUNVLS)
!             Input parameters:
!               DIMENS = 1
!               X(1)      The coordinate of the evaluation point.
!               NUMFUN Integer that defines the number of
!                      components of I.
!             Output parameter:
!               FUNVLS Real array of dimension NUMFUN
!                      that defines NUMFUN components of the integrand.
!   VER    Real array of dimension (1,2).
!          The coordinates of the vertices of the interval.
!   KEY    Integer
!          key for choice of local integration rule
!          a gauss-kronrod pair is used with
!           7 - 15 points if key < 2,
!          10 - 21 points if key = 2,
!          15 - 31 points if key = 3,
!          20 - 41 points if key = 4,
!          25 - 51 points if key = 5,
!          30 - 61 points if key > 5.
!
!   ON RETURN
!
!   BASVAL Real array of dimension NUMFUN.
!          The values for the basic rule for each component
!          of the integrand.
!   RGNERR Real array of dimension NUMFUN.
!          The error estimates for each component of the integrand.
!
!***references  (none)
!***routines called  Integrand
!***end prologue  dqknn
! 
      REAL(kind=stnd), DIMENSION(:,:), INTENT(IN) :: VER
      REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: wg,xgk,wgk
      INTEGER, INTENT(IN) :: NUMFUN,KEY
      REAL(kind=stnd), DIMENSION(:), INTENT(OUT) :: BASVAL,RGNERR
      INTERFACE 
         FUNCTION Integrand(NUMFUN,X) RESULT(Value)
            USE Precision_Model
            INTEGER, INTENT(IN) :: NUMFUN
            REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
            REAL(kind=stnd), DIMENSION(NUMFUN) :: Value
         END FUNCTION Integrand
      END INTERFACE

      REAL(kind=stnd):: absc,centr,dhlgth,hlgth,resasc, reskh
      REAL(kind=stnd), DIMENSION(NUMFUN) :: fc, fsum,fval1,fval2, &
                                            resabs, resg, resk
      REAL(kind=stnd), DIMENSION(size(xgk)-1,NUMFUN) ::  fv1,fv2
      INTEGER ::  j,jtw,jtwm1,m,upper,dummy
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 15-point kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 7-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 7-point gauss rule
!
!           wgk    - weights of the 15-point kronrod rule
!
!           wg     - weights of the 7-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.
!
!           list of major variables
!           -----------------------
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 7-point gauss formula
!           resk   - result of the 15-point kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
!***first executable statement  dqk15
!
   centr = 0.5_stnd*(ver(1,1)+ver(1,2))
   hlgth = 0.5_stnd*(ver(1,2)-ver(1,1))
   dhlgth = abs(hlgth)
!
!           compute the 15-point kronrod approximation to
!           the integral, and estimate the absolute error.
!
   upper = size(wgk)
   fc = Integrand(NUMFUN,(/centr/))
   DO m=1,NUMFUN
      SELECT CASE (key)
      CASE (:1,3,5)
         resg(m) = fc(m)*wg(size(wg))
      CASE (2,4,6:)
         resg(m) = 0.0_stnd
      END SELECT
      resk(m) = fc(m)*wgk(upper)
      resabs(m) = abs(resk(m))
   END DO
   dummy = (upper - 1)/2
   DO j=1,dummy
      jtw = j*2
      absc = hlgth*xgk(jtw)
      fval1 = Integrand(NUMFUN,(/centr-absc/))
      fval2 = Integrand(NUMFUN,(/centr+absc/))
      fv1(jtw,:) = fval1
      fv2(jtw,:) = fval2
      fsum = fval1+fval2
      DO m=1,NUMFUN
         resg(m) = resg(m)+wg(j)*fsum(m)
         resk(m) = resk(m)+wgk(jtw)*fsum(m)
         resabs(m) = resabs(m)+wgk(jtw)*(abs(fval1(m))+abs(fval2(m)))
      END DO
   END DO
   dummy = size(wg)
   DO j = 1,dummy
      jtwm1 = j*2-1
      absc = hlgth*xgk(jtwm1)
      fval1 = Integrand(NUMFUN,(/centr-absc/))
      fval2 = Integrand(NUMFUN,(/centr+absc/))
      fv1(jtwm1,:) = fval1
      fv2(jtwm1,:) = fval2
      fsum = fval1+fval2
      DO m=1,NUMFUN
         resk(m) = resk(m)+wgk(jtwm1)*fsum(m)
         resabs(m) = resabs(m)+wgk(jtwm1)*(abs(fval1(m))+abs(fval2(m)))
      END DO
   END DO
   DO m=1,NUMFUN
      reskh = resk(m)*0.5_stnd
      resasc = wgk(upper)*abs(fc(m)-reskh)
      DO j=1,upper-1
        resasc = resasc+wgk(j)*(abs(fv1(j,m)-reskh)+abs(fv2(j,m)-reskh))
      END DO
      BASVAL(m) = resk(m)*hlgth
      resabs(m) = resabs(m)*dhlgth
      resasc = resasc*dhlgth
      RGNERR(m) = abs((resk(m)-resg(m))*hlgth)
      IF (resasc /= 0.0_stnd .AND. RGNERR(m) /= 0.0_stnd) THEN
           RGNERR(m) = resasc*min(1.0_stnd,(200.0_stnd*RGNERR(m)/resasc)**1.5_stnd)
      END IF
      IF (resabs(m) > TINY(RGNERR(m))/(50.0_stnd*EPSILON(RGNERR(m))))&
         THEN
            RGNERR(m) = max((EPSILON(RGNERR(m))*50.0_stnd)*resabs(m),RGNERR(m))
      END IF
   END DO
   RETURN
END SUBROUTINE Dqknn

END Module QuadratureRule
