function kda_top () {
var path=kpath + 'aps/';
var ts=new Date().getTime();
var ms=new Date().getMilliseconds();
var Top_nam=new Array(); var Top_ext=new Array(); var Top_wid=new Array(); var Top_hgt=new Array();
var Top_alt=new Array(); var Top_txt=new Array(); var Top_url=new Array(); var Top_wgt=new Array(); var wgt=0;
wgt+=1; Top_wgt[1]=wgt;
Top_nam[1]="t-drb-19m04-5ai"; Top_ext[1]=".jpg"; Top_wid[1]=750; Top_hgt[1]=100;
Top_url[1]='https://www.datarobot.com/lp/risk/?utm_source=kdnuggets&utm_medium=banner&utm_campaign=5SOLUTIONSRISKkdn';
Top_txt[1]="5 AI Solutions Every Chief Risk Officer Needs - See what you are missing";
Top_alt[1]="txt";
wgt+=1; Top_wgt[2]=wgt;
Top_nam[2]="t-jmp-19m01-berry"; Top_ext[2]=".jpg"; Top_wid[2]=750; Top_hgt[2]=100;
Top_url[2]='https://www.jmp.com/en_us/offers/data-mining-techniques-book.html?utm_source=kdnuggetshp&utm_medium=advertisement&utm_campaign=wp701a0000000t5sD';
Top_txt[2]="<font size=+1>Derived Variables: Making the data mean more - Download a free book chapter</font>";
Top_alt[2]="JMP";
wgt+=0.5; Top_wgt[3]=wgt;
Top_nam[3]="t-paw-19m04-dlwde"; Top_ext[3]=".jpg"; Top_wid[3]=750; Top_hgt[3]=100;
Top_url[3]='https://deeplearningworld.de/?utm_source=kdnbanner';
Top_txt[3]="Deep Learning World: The premier conference on deployment of Deep Learning. Munich, 6-7 May 2019";
Top_alt[3]="txt";
wgt+=0.5; Top_wgt[4]=wgt;
Top_nam[4]="t-paw-19m04-indde"; Top_ext[4]=".jpg"; Top_wid[4]=750; Top_hgt[4]=100;
Top_url[4]='https://predictiveanalyticsworld.de/en/industry4-0/muenchen2019/?utm_source=kdnbanner';
Top_txt[4]="PAW Industry 4.0: Putting <b>Machine Intelligence</b> Into Production. Munich, 6-7 May 2019";
Top_alt[4]="txt";
wgt+=1; Top_wgt[5]=wgt;
Top_nam[5]="t-iot-19m04-duc"; Top_ext[5]=".png"; Top_wid[5]=750; Top_hgt[5]=100;
Top_url[5]='https://io-tahoe.com/';
Top_txt[5]="Discover, understand, and catalog your data with Io-Tahoe. Smart Data Discovery - AI-driven Catalog";
Top_alt[5]="txt";
wgt+=1; Top_wgt[6]=wgt;
Top_nam[6]="t-car-19m04"; Top_ext[6]=".gif"; Top_wid[6]=1; Top_hgt[6]=1;
Top_url[6]='<ins class="dcmads" style="display:inline-block;width:750px;height:100px" data-dcm-placement="N5506.2008108KDNUGGETS/B22398349.242178943" data-dcm-rendering-mode="iframe" data-dcm-https-only data-dcm-resettable-device-id="" data-dcm-app-id=""> <script src="https://www.googletagservices.com/dcm/dcmads.js"></script></ins>';
Top_txt[6]="Intuit Careers - Learn More";
Top_alt[6]="Intuit";
var rtop=Math.random()*5;
var n_ad=1; while (Top_wgt[n_ad]<rtop) {n_ad++};
var out='<table border=0 cellspacing=5 cellpadding=20 width=990><tr><td valign=top align=center width=750>';
  if (Top_wid[n_ad] == 1) {
      out+=Top_url[n_ad];
  } else {
      out+='<a href=' + Top_url[n_ad] + ' onclick="ga(\'send\',\'pageview\',\'/zt/' + Top_nam[n_ad] + '\');" target=_blank>';
  }
  var whb=' width=' + Top_wid[n_ad] + ' height=' + Top_hgt[n_ad] + ' border=0';
  out+='<img src=' + path + Top_nam[n_ad] + Top_ext[n_ad] + '?ms=' + ms + whb + ' alt="' + Top_alt[n_ad] + '">';
  out+='<br><b>' + Top_txt[n_ad] + '</b></a></td>';
  out+='</tr></table><br>';document.writeln(out);}
function kda_sid_init() {
var ts=new Date().getTime();
Sid_nam=new Array(); Sid_ext=new Array(); Sid_txt=new Array(); Sid_url=new Array();
Sid_alt=new Array(); Sid_wid=new Array(); Sid_hgt=new Array(); Sid_wgt=new Array(); var wgt=0;
wgt+=0.89; Sid_wgt[1]=wgt;
Sid_nam[1]="e-sas-19m01-an"; Sid_ext[1]=".gif"; Sid_wid[1]=1; Sid_hgt[1]=1;
Sid_url[1]='<ins class="dcmads" style="display:inline-block;width:300px;height:250px" data-dcm-placement="N6626.289580.KDNUGGETS.COM/B22085569.236894257" data-dcm-rendering-mode="iframe" data-dcm-https-only data-dcm-resettable-device-id="" data-dcm-app-id=""><script src="https://www.googletagservices.com/dcm/dcmads.js"></script></ins>';
Sid_txt[1]=" ";
Sid_alt[1]="SAS Analytics Campaign";
wgt+=0.89; Sid_wgt[2]=wgt;
Sid_nam[2]="e-que-19m02-mma"; Sid_ext[2]=".jpg"; Sid_wid[2]=300; Sid_hgt[2]=250;
Sid_url[2]='https://smith.queensu.ca/grad_studies/mma/landing.php?utm_source=KDNuggets&utm_medium=email&utm_campaign=GMMA2020';
Sid_txt[2]="Smith School<br>Master of Management Analytics<br> Unleash the True Potential of Data";
Sid_alt[2]="Smith School Master of Management Analytics<br> Unleash the True Potential of Data";
wgt+=0.89; Sid_wgt[3]=wgt;
Sid_nam[3]="e-dom-19m01-rev2"; Sid_ext[3]=".png"; Sid_wid[3]=300; Sid_hgt[3]=250;
Sid_url[3]='https://rev.dominodatalab.com/?utm_source=kdnuggets&utm_medium=display&utm_campaign=KDNuggets%20Paid_CH';
Sid_txt[3]="REV 2 Data Science Leaders Summit, May 23-24, New York City<br>Get $100 off your order by using promo code: KDNuggetREV";
Sid_alt[3]="REV 2 Data Science Leaders Summit, May 23-24, NYC";
wgt+=0.44; Sid_wgt[4]=wgt;
Sid_nam[4]="e-nyu-18m08-b"; Sid_ext[4]=".gif"; Sid_wid[4]=300; Sid_hgt[4]=250;
Sid_url[4]='https://web-marketing.stern.nyu.edu/global-programs/business-analytics/?utm_source=kdnuggets&utm_medium=ros&utm_campaign=kdn12919&utm_term=300x250&utm_content=fy19';
Sid_txt[4]="How you view data changes how you view business decisions. Take the next step. Get MSBA at NYU Stern.";
Sid_alt[4]="NYU";
wgt+=0.44; Sid_wgt[5]=wgt;
Sid_nam[5]="e-nyu-18m08-c"; Sid_ext[5]=".gif"; Sid_wid[5]=300; Sid_hgt[5]=250;
Sid_url[5]='https://web-marketing.stern.nyu.edu/global-programs/business-analytics/?utm_source=kdnuggets&utm_medium=ros&utm_campaign=kdn12919&utm_term=300x250&utm_content=fy19';
Sid_txt[5]="How you look at data changes how you look at business strategies. Take the next step. Get MSBA at NYU Stern.";
Sid_alt[5]="NYU";
wgt+=0.89; Sid_wgt[6]=wgt;
Sid_nam[6]="e-gur-19m04"; Sid_ext[6]=".jpg"; Sid_wid[6]=300; Sid_hgt[6]=250;
Sid_url[6]='http://www.gurobi.com/products/gurobi-optimizer';
Sid_txt[6]="The Fastest<br>Mathematical<br>Optimization Solver<br>Try It Now!";
Sid_alt[6]="The Fastest Mathematical<br>Optimization Solver<br>Try It Now!";
wgt+=0.89; Sid_wgt[7]=wgt;
Sid_nam[7]="e-car-19m04-intu"; Sid_ext[7]=".gif"; Sid_wid[7]=1; Sid_hgt[7]=1;
Sid_url[7]='<ins class="dcmads" style="display:inline-block;width:300px;height:250px" data-dcm-placement="N5506.2008108KDNUGGETS/B22398349.242073350" data-dcm-rendering-mode="iframe" data-dcm-https-only data-dcm-resettable-device-id="" data-dcm-app-id=""> <script src="https://www.googletagservices.com/dcm/dcmads.js"></script></ins>';
Sid_txt[7]="Intuit Careers - Learn More";
Sid_alt[7]="Intuit";
wgt+=0.89; Sid_wgt[8]=wgt;
Sid_nam[8]="e-kni-19m04-fall"; Sid_ext[8]=".jpg"; Sid_wid[8]=300; Sid_hgt[8]=250;
Sid_url[8]='https://www.knime.com/about/events/knime-fall-summit-2019-austin?kd0319';
Sid_txt[8]="KNIME Fall Summit 2019<br>Nov 5-8, Austin<br>Use code KDNUGGETS for 10% off";
Sid_alt[8]="text";
wgt+=0.44; Sid_wgt[9]=wgt;
Sid_nam[9]="e-tdw-19m04-ch-ds-reg"; Sid_ext[9]=".jpg"; Sid_wid[9]=300; Sid_hgt[9]=250;
Sid_url[9]='https://tdwi.org/events/conferences/chicago/home.aspx?utm_source=KD&utm_medium=banner&utm_campaign=chicago';
Sid_txt[9]="TDWI Chicago, Apr 28 - May 3<br>Full Day Data Strategy Course<br>Save your seat";
Sid_alt[9]="text";
wgt+=0.44; Sid_wgt[10]=wgt;
Sid_nam[10]="e-tdw-19m04-ch-ac-reg"; Sid_ext[10]=".jpg"; Sid_wid[10]=300; Sid_hgt[10]=250;
Sid_url[10]='https://tdwi.org/events/conferences/chicago/home.aspx?utm_source=KD&utm_medium=banner&utm_campaign=chicago';
Sid_txt[10]="TDWI Chicago, Apr 28 - May 3<br>Tame Your Analytics Chaos<br>Save your seat";
Sid_alt[10]="text";
wgt+=0.9; Sid_wgt[11]=wgt;
Sid_nam[11]="e-psu-19m03-ubd"; Sid_ext[11]=".jpg"; Sid_wid[11]=300; Sid_hgt[11]=250;
Sid_url[11]='https://www.worldcampus.psu.edu/degrees-and-certificates/business-analytics-certificate/overview?utm_source=datasciencecentral&utm_medium=banner&utm_campaign=DATA-BA+18-19&utm_content=300x250&cid=BNNR35117';
Sid_txt[11]="Earn a Business Analytics<br>Certificate online<br>at Penn State";
Sid_alt[11]="text";
}
function kda_sid_write (nads) {
  var path=kpath + 'aps/';
  var ms=new Date().getMilliseconds(); var adshown=new Array();
  for (adpos=1; adpos<=nads; adpos++) {
    do { var adn=1; var re=Math.random()*8;
      while (Sid_wgt[adn]<re) {adn++};
    } while (adshown[adn]==1);
    if (Sid_wid[adn]==1) { s=Sid_url[adn]; s=s.replace('ze/e-','ze/' + nads + adpos + ':e-'); 
      s+='<img src=' + path + Sid_nam[adn] + '.gif?ms' + ms + ' width=1 height=1 border=0>';
    } else {
      s='<a href=' + Sid_url[adn] + ' onclick="ga(\'send\',\'pageview\',\'/ze/' + nads + adpos + ':' + Sid_nam[adn] + '\');" target=_blank>';
      s+='<img src=' + path + Sid_nam[adn] + Sid_ext[adn] + '?ms' + ms +  ' width=' + Sid_wid[adn] + ' height=' + Sid_hgt[adn] + ' border=0 alt="' + Sid_alt[adn] + '">';
      s+='<br><b>' + Sid_txt[adn] + '</b></a><br>';}
    if (nads>1) {s+='<br><hr class="grey-line"><br>'};
    document.writeln(s); adshown[adn]=1;
  }
}
