function $id(id) { return document.getElementById(id); }
var UI = {
  chkmsgfr:'Verifiez que la calculatrice Numworks est connectee.\nCliquez sur le bouton Detecter puis reessayez.',
  chkmsgen:'Check that the Numworks calculator is connected.\nClick on the detect button and try again.',
  xcasfrmsg:"Les liens Test ci-dessus chargent vos scripts dans MicroPython de Xcas. Vous pouvez ensuite les &eacute;diter, les envoyer par email, les poster sur des forums (par exemple <a href=\"https://xcas.univ-grenoble-alpes.fr/forum/\">le forum de Xcas</a>) ou les renvoyer vers la Numworks (icones enveloppe, boutons F et Numworks→ en haut). Pour visualiser un graphique Kandinsky, matplotlib ou tortue validez ; ou , ou . dans une ligne de commande vide.",
  xcasenmsg:"The Test link above will load your scripts into Xcas's Micropython. From there you can edit, execute, send them by mail, post them on a forum (for eaemple <a href=\"https://xcas.univ-grenoble-alpes.fr/forum/\">Xcas forum</a>) and send them back to the Numworks (from the toolbar at the top, mail icon, F and  Numworks→ buttons). If you want to see a Kandinsky, matplotlib or turtle graph, enter ; or , or . on an empty commandline and hit Enter",
  nws_progress:0,
  nws_progresslegend:0,
  Datestart:Date.now(),
  langue: -1, // -1 french
  calc: 2, // 1 KhiCAS, 2 Numworks, 3 TI Nspire CX
  calculator:0, // !=0 if hardware Numworks connected
  calculator_connected:false,
  tar_first_modif_offset:0,
  nws_connect:function(){
    UI.nws_progress=$id('plog');
    UI.nws_progresslegend=$id('ploglegend');
    if (navigator.usb){
      //console.log('nws_connect 0');
      UI.calculator=0;
      UI.calculator_connected=false;
      //console.log('nws_connect 1');
      function autoConnectHandler(e) {
	UI.calculator.stopAutoConnect();
	console.log('connected');
	UI.calculator_connected=true;
      }
      UI.calculator= new Numworks();
      console.log('nws_connect',UI.calculator);
      navigator.usb.addEventListener("disconnect", function(e) {
	if (UI.calculator==0) return;
	UI.calculator.onUnexpectedDisconnect(e, function() {
	  UI.calculator_connected=false;
          UI.calculator=0;
	});
      });
      //console.log('nws_connect 2');
      UI.calculator.autoConnect(autoConnectHandler);
      return UI.calculator;
      //console.log('nws_connect 3');
    }
  },
  nws_detect_success:function(){
    alert('Success');
  },
  nws_detect_failure:function(error){
    console.log(error);
    alert(UI.langue==-1?'Verifiez la connexion de la calculatrice et reessayez':'Check calculator connection and try again');
  },
  nws_detect:async function(success,failure){
    UI.calculator= new Numworks();
    UI.calculator.detect(success,failure);
  },
  nws_rescue_connect:async function(){
    if (navigator.usb){
      //console.log('nws_rescue_connect 0');
      UI.calculator=0;
      UI.calculator_connected=false;      //console.log('nws_rescue_connect 1');
      function autoConnectHandler(e) {
	UI.calculator.stopAutoConnect();
	console.log('rescue connected');
	UI.calculator_connected=true;
      }
      UI.calculator= new Numworks.Recovery();
      console.log('nws_rescue_connect',UI.calculator);
      navigator.usb.addEventListener("rescue disconnect", function(e) {
	if (UI.calculator==0) return;
	UI.calculator.onUnexpectedDisconnect(e, function() {
	  UI.calculator_connected=false;
          UI.calculator=0;
	});
      });
      //console.log('nws_rescue_connect 2');
      if (await UI.calculator.autoConnect(autoConnectHandler)==-1){
	return -1;
      }
      //console.log('nws_rescue_connect 3');
      return 0;
    }
  },
  adjust_exam_sector:function(data){
    let l=data.byteLength;
    if (l<0x3000) return;
    let tab=new Uint8Array(data);
    for (let i=0x1000;i<0x2000;i++)
      tab[i]=0xff;
    for (let i=0x2000;i<0x3000;i++)
      tab[i]=0x0;
  },
  sig_check:async function(sig,data,fname){
    // sig should be a list of lists of size 3 (name, length, hash)
    let i=0,l=sig.length;
    for (;i<l;++i){
      let cur=sig[i];
      if (cur[0]!=fname) continue;
      console.log('sig_check',i,cur[1],data.byteLength);
      if (cur[1]>data.byteLength) continue;
      let dat=data.slice(0,cur[1]);
      let digest = await window.crypto.subtle.digest('SHA-256', dat);
      digest=Array.from(new Uint8Array(digest));
      console.log(cur[2],digest);
      let j=0;
      for (;j<32;++j){
	let tst=(digest[j]-cur[2][j]) % 256;
	// console.log(j,digest[j],cur[2][j]);
	if (tst)
	  break;
      }
      if (j==32){
	console.log('signature match',cur[0]);
	return true;
      }      
    }
    return false;
  },
  numworks_load:function(backup=false){
    UI.calc=2;
    UI.nws_connect();
    window.setTimeout(UI.numworks_load_,100,backup);
  },
  backup_load:function(file){
    let reader = new FileReader();
    reader.onerror = function (evt) { }
    reader.readAsArrayBuffer(file);
    reader.onload = function (evt) {
      //UI.backup2rec(evt.target.result);
      UI.storage_records=UI.backup2rec(evt.target.result);
      UI.store2html(UI.storage_records,'storelist');
    }
  },
  backup2rec:function(buf){
    let pos=4,nwstoresize=32768;
    let view=new Uint8Array(buf);
    let rec=[];
    for (;pos<nwstoresize;){
      let L=view[pos+1]*256+view[pos];
      pos += 2;
      if (L==0) break;
      L -= 2;
      let name='',type='',point=false,i;
      for (i=0;i<L;++i){
	let ch=view[pos+i];
	if (ch==0) break;
	if (String.fromCharCode(ch)=='.'){
	  point=true;
	  continue;
	}
	if (point)
	  type += String.fromCharCode(ch);
	else
	  name += String.fromCharCode(ch);
      }
      let len=L-i-2;
      let data=new Uint8Array(buf,pos+i+1,len);
      let code="";
      if (type=="py")
	code=new TextDecoder().decode(new Uint8Array(buf,pos+i+2,len-1));
      let b=new Blob([data], {type: 'application/octet-binary'});
      let record={"name":name, "type":type,  "autoImport": false, "data": b,"code":code,"size": len}
      //console.log(record);
      rec.push(record)
      pos += L;
    }
    return rec;
  },
  numworks_retry:function(){
    alert(UI.langue==-1?'Calculatrice detectee. Reessayez!':'Calculator detected! Try again');
  },
  numworks_load_: async function(backup){
    console.log(UI.calculator,UI.calculator_connected);
    if (UI.calculator==0 || !UI.calculator_connected){
      if (UI.calculator) UI.calculator.stopAutoConnect();
      UI.nws_detect(UI.numworks_retry,UI.nws_detect_failure);
      return;
    }
    if (backup){
      let pinfo = await UI.calculator.getPlatformInfo();
        
      let storage_blob = await UI.calculator.__retreiveStorage(pinfo["storage"]["address"], pinfo["storage"]["size"]);
      filename = "backup.nws";
      saveAs(storage_blob.slice(0,32768), filename);
      return;
    }
    let storage = await UI.calculator.backupStorage();
    let rec=UI.storage_records=storage.records;
    let ht=UI.store2html(rec,'storelist'); //console.log(ht);
    return;
    let j=0,s='Choisissez un numero parmi ';
    for (;j<rec.length;++j){
      s+=j;
      s+=':'+rec[j].name+', ';
    }
    let p=prompt(s),n=0;
    if (!p) return;
    let l=p.length;
    for (j=0;j<l;++j){
      if (p[j]<'0' || p[j]>'9'){ alert(UI.langue==-1?'Nombre invalide':'Invalid number'); return ;}
      n*=10;
      n+=p.charCodeAt(j)-48;
    }
    if (n>=rec.length){ alert(UI.langue==-1?'Choix invalide':'Invalide choice'); return; }
    s=rec[n].code; p=rec[n].name;
    let blob = new Blob([s], {type: "text/plain;charset=utf-8"});
    filename = p + ".py";
    saveAs(blob, filename);
  },
  loadFile:function(file){
    let reader = new FileReader();
    reader.onerror = function (evt) { }
    let backup=file.name.length>4 && file.name.substr(file.name.length-4,4)==".nws";
    if (backup){
      if (!confirm(UI.langue==-1?'Remplacer tous les scripts de la calculatrice?':'Overwrite all calculator scripts?'))
	return;
      reader.readAsArrayBuffer(file);
    }
    else
      reader.readAsText(file, "UTF-8");
    reader.onload = function (evt) {
      let data=evt.target.result;
      if (!backup)
	data+=String.fromCharCode(0);
      UI.numworks_save(file.name,data);
    }
  },
  numworks_save:function(filename,S){
    //console.log(filename,S); return;
    UI.nws_connect();
    //window.setTimeout(UI.numworks_save_,100,filename,S);
    window.setTimeout(UI.numworks_save_,100,filename,S);
  },
  numworks_save_:async function(filename,S){
    if (UI.calculator==0 || !UI.calculator_connected){
      alert(UI.langue==-1?UI.chkmsgfr:chkmsgen);
      if (UI.calculator) UI.calculator.stopAutoConnect();
      return -1;
    }
    if (filename.length>4 && filename.substr(filename.length-4,4)==".nws"){
      let pinfo = await UI.calculator.getPlatformInfo();
      UI.calculator.device.startAddress = pinfo["storage"]["address"];
      //console.log(UI.calculator.device.startAddress,S); return;
      let res=await UI.calculator.device.do_download(UI.calculator.transferSize, S, false);
      if (res==0)
	alert('Unable to store to the calculator');
      return res;
    }
    if (UI.is_py(filename))
      filename=filename.substr(0,filename.length-3);
    let storage = await UI.calculator.backupStorage();
    let rec=storage.records,j=0;
    for (;j<rec.length;++j){
      if (rec[j].name==filename){
	if (!confirm((UI.langue==-1?'? Ecraser ':'? Overwrite ')+rec[j].name))
	  return;
	rec[j].code=S;
	break;
      }
    }
    if (j==rec.length)
      rec.push({"name": filename, "type":"py", "autoImport": false, "code": S});
    await UI.calculator.installStorage(storage, function() {
      // Do stuff after writing to the storage is done
      console.log(filename+'.py saved to Numworks');
    });
    return 0;
  },    
  numworks_install_delta:function(do_backup=1,khionly=0){
    if (khionly && !confirm(UI.langue==-1?'En confirmant l\'installation, vous acceptez la licence CC BY-NC-SA de Khi':'If you confirm, you accept the CC BY-NC-SA license of Khi.'))
      return;
    if (!khionly && !confirm(UI.langue==-1?'En confirmant l\'installation, vous acceptez les licences CC BY-NC-SA de Khi et GPL2 de KhiCAS, BSD2 Hexedit, LGPL2 de Nofrendo, MIT de Peanut-GB':'If you confirm install, you accept the CC BY-NC-SA license of Khi and GPL2 license of KhiCAS, BSD2 for Hexedit, LGPL2 for Nofrendi, MIT for Peanut-GB'))
      return;
    $id('persoinit').style.display='none';    
    UI.calc=2;
    UI.nws_connect();
    window.setTimeout(UI.numworks_install_delta_,100,do_backup,khionly);
  },
  numworks_install_delta_:async function(do_backup,khionly){
    if (khionly==5){
      UI.calculator.device.startAddress = 0x90180000;
      UI.calculator.device.logProgress(0,'khi slot B');
      let data=await UI.loadfile('khi.B.bin');
      let res=await UI.calculator.device.do_download(UI.calculator.transferSize, data, false);
      if (res==0){
	alert(UI.langue==-1?'Erreur d\'ecriture flash. Tapez reset+4 et reessayez':'Error writing flash. Type reset+4 and try again.');
	return 3;
      }
      alert(UI.langue==-1?'Tapez reset+2 pour booter un lanceur minimal de KhiCAS sur le slot B. Ceci augmente la memoire de 50% environ. Attention, cela ne fonctionnera que si vous installez le bootloader de cette page!':'Type reset+2 to boot a minimal KhiCAS launcher from slot B. This increases memory available. Beware this works only if you install the bootloader of this page!');
      UI.print('Khi slot B OK, press RESET on calculator back');
      return 0;
    }
    if (khionly==4){
      UI.calculator.device.startAddress = 0x08000000;
      UI.calculator.device.logProgress(0,'protection');
      let data=await UI.loadfile('bootloader.bin');
      let res=0;
      res=await UI.calculator.device.do_download(UI.calculator.transferSize, data, false);
      if (res==0){
	alert(UI.langue==-1?'Erreur d\'ecriture flash interne. Tapez reset+6 et reessayez. Si reset+6 ne fonctionne pas, essayez de deverrouiller votre Numworks sur phi.getomega.dev':'Error writing internal flash. Type reset+6 and try again. If reset+6 does not work, try to unlock your Numworks on phi.getomega.dev');
	return 3;
      }
      alert(UI.lang==-1?'Le multi-boot permet de lancer 2 firmwares A (<=1.5M) et B(<=0.5M) : reset-1 ou reset-2 permet de changer. Si vous installez Khi A et B, votre calculatrice sera protegee contre des mises a jour (risque de verrouillage). Cette protection se desactive temporairement par reset-4, par exemple pour mettre a jour Khi. Tapez reset+6 pour effacer le multiboot':'With multi-boot installed you can run 2 firmwares A (<=1.5M) and B (<=0.5M). If you install Khi A and B, your calculator will be protected against upgrade (lock risk). Type reset-4 to temporarily remove firmware protection, e.g. to upgrade Khi. Type reset-6 to remove the multi-boot');
      UI.print('Bootloader OK, press RESET on calculator back');
      return 0;
    }
    if (UI.calculator==0 || !UI.calculator_connected){
      alert(UI.langue==-1?UI.chkmsgfr:chkmsgen);
      if (UI.calculator) UI.calculator.stopAutoConnect();
      return -1;
    } 
    UI.print(khionly?'=========== installing Khi':'=========== installing Khi and KhiCAS');
    let data,res;
    if (!khionly){
      UI.print('erase/write KhiCAS apps, wait about 2 minutes');
      data=await UI.loadfile('apps.tar');
      UI.calculator.device.startAddress = 0x90200000;
      UI.calculator.device.logProgress(0,'apps');
      res=await UI.calculator.device.do_download(UI.calculator.transferSize, data, false);
      UI.print(res?'KhiCAS installed':'error installing KhiCAS');
      if (res==0){
	alert('error installing KhiCAS');
	return 1;
      }
    }
    let installboot=0;
    let pinfo = await UI.calculator.getPlatformInfo();
    if (pinfo["bootloader"]==0 || khionly==3){
      installboot=1;
      if (khionly==3) // skip protection question parameter
	khionly=1;
      else {
	if (1 ||
	    confirm(UI.langue==-1?'Votre calculatrice n\'est pas protegee contre un verrouillage. Voulez-vous la proteger?':'Your calculator is not protected against an hostile lock. Protect?')){
	  alert(UI.langue==-1?'Le bootloader de protection sera installe. Vous pourrez l\'effacer avec reset+6. Attention, ceci ne protege pas d\'une mise a jour sur le site de Numworks.':'Protecting bootloader will be installed, you will be able to remove it with reset+6. Beware, this does not protect from running an upgrade from Numworks.');
	  khionly=2;
	}
	else
	  khionly=1;
      }
    }
    else
      khionly=2;
    UI.print('erase/write Khi external, wait about 1/2 minute');
    UI.calculator.device.startAddress = 0x90000000;
    let firmwarename=khionly==2?'khi.A.bin':'delta.external.bin';
    data=await UI.loadfile(firmwarename);
    UI.calculator.device.logProgress(0,firmwarename);
    res=await UI.calculator.device.do_download(UI.calculator.transferSize, data, false);
    if (res==0){
      alert('Error writing Khi');
      return 2;
    }
    UI.calculator.device.startAddress = 0x90180000;
    UI.calculator.device.logProgress(0,'khi slot 2');
    data=await UI.loadfile('khi.B.bin');
    res=await UI.calculator.device.do_download(UI.calculator.transferSize, data, false);
    if (installboot){
      UI.print('Khi external OK, erase/write internal');
      UI.calculator.device.startAddress = 0x08000000;
      UI.calculator.device.logProgress(0,'internal');
      data=await UI.loadfile(khionly==2?'bootloader.bin':'delta.internal.bin');
      res=await UI.calculator.device.do_download(UI.calculator.transferSize, data, false);
      if (res==0){
	alert(UI.langue==-1?'Erreur d\'ecriture flash interne. Tapez reset+6, installez le mode de recuperation et reessayez.':'Error writing internal flash. Type reset+6, run rescue mode and try again');
	return 3;
      }
      UI.print((installboot?'Bootloader OK':'Khi internal OK')+', press RESET on calculator back');
    }
    else 
      UI.print('Khi external OK');
    alert(UI.langue==-1?'Installation terminee. Appuyer sur RESET a l\'arriere de la calculatrice':'Install success. Press RESET on the calculator back.');
    return 0;
  },
  loadfile:function(file) {
    return fetch(file).then(function(response) {
      if (!response.ok) {
        throw new Error("HTTP error, status = " + response.status);
      }
      return response.arrayBuffer();
    })
  },
  numworks_rescue:async function(){
    alert(UI.langue==-1?'Connectez la calculatrice, appuyez sur la touche 6 de la calculatrice, enfoncez un stylo dans le bouton RESET au dos en laissant la touche 6 appuyee, relachez la touche 6, l\'ecran doit etre eteint et la diode rouge allumee':'Connect the calculator, press the 6 key on the calculator, press the RESET button on the back keeping the 6 key pressed, release the 6 key, the screen should be down and the led should be red');
    UI.calc=2;
    let res=await UI.nws_rescue_connect();
    if (res==-1){
      return;
    }
    window.setTimeout(UI.numworks_rescue_,1000);
  },
  numworks_rescue_:async function(sigfile,rwcheck){
    if (UI.calculator==0 || !UI.calculator_connected){
      alert(UI.langue==-1?UI.chkmsgfr:chkmsgen);
      //console.log('numworks_rescue_',UI.calculator);
      if (UI.calculator) UI.calculator.stopAutoConnect();
      UI.calculator=0;
      return -1;
    }
    UI.print(UI.langue==-1?'Envoi du mode de recuperation, patientez environ 15 secondes':'Sending rescue mode, wait about 15 secondes');
    UI.calculator.device.startAddress = 0x20030000;
    let data=await UI.loadfile('recovery');
    await UI.calculator.device.clearStatus();
    let res=await UI.calculator.device.do_download(UI.calculator.transferSize, data, true);
    if (res==0){
      alert('Error sending recovery');
      return 1;
    }
    let tmp=confirm(UI.langue==-1?'La calculatrice devrait etre en mode recuperation. Cliquer sur le bouton Detecter puis sur install KhiCAS':'Calculator should be in rescue mode. Click on the Detect button then install KhiCAS.');
    return 0;
    if (tmp)
      UI.nws_detect(UI.numworks_install_delta,UI.nws_detect_failure);
    else
      UI.nws_detect(UI.numworks_detect_success,UI.nws_detect_failure);
  },
  numworks_buffer:0,
  read_string:function(buffer,str_offset, size) {
    let strView = new Uint8Array(buffer, str_offset, size);
    let i = 0;
    let rtnStr = "";
    while(strView[i] != 0) {
      rtnStr += String.fromCharCode(strView[i]);
      i++;
    }
    return rtnStr;
  },
  // tar support adapted from tarball.js https://github.com/ankitrohatgi/tarballjs (MIT licence)
  filename(buffer,header_offset) {
    let name = UI.read_string(buffer,header_offset, 100);
    return name;
  },
  filetype(buffer,header_offset) {
    // offset: 156
    let typeView = new Uint8Array(buffer, header_offset+156, 1);
    let typeStr = String.fromCharCode(typeView[0]);
    if(typeStr == "0") {
      return "file";
    } else if(typeStr == "5") {
      return "directory";
    } else {
      return typeStr;
    }
  },
  filesize(buffer,header_offset) {
    // offset: 124
    let szView = new Uint8Array(buffer, header_offset+124, 12);
    let szStr = "";
    for(let i = 0; i < 11; i++) {
      let tmp=szView[i];
      if (tmp<48 || tmp>57) return -1; // invalid file size
      szStr += String.fromCharCode(tmp);
    }
    return parseInt(szStr,8);
  },
  numworks_chk_buffer:function(buffer){
    if (buffer==0){
      alert(UI.langue==-1?'Commencez par recuperer la flash de la calculatrice':'First get the calculator flash');
      return 0;
    }
    // FIXME should check that buffer is valid
    return 1;
  },
  numworks_tarinfo:function(buffer){
    if (!UI.numworks_chk_buffer(buffer))
      return [];
    let offset = 0;
    let file_size = 0;       
    let file_name = "";
    let file_type = null;
    let file_info=[];
    while(offset < buffer.byteLength - 512) {
      file_name = UI.filename(buffer,offset); // file name
      if (file_name.length == 0) break;
      file_type = UI.filetype(buffer,offset);
      file_size = UI.filesize(buffer,offset);
      if (file_size<0) break;
      file_info.push({
        "name": file_name,
        "type": file_type,
        "size": file_size,
        "header_offset": offset
      });
      offset += (512 + 512*Math.trunc(file_size/512));
      if(file_size % 512) {
        offset += 512;
      }
    }
    return file_info;
  },
  leftpad:function(number,targetLength) {
    let output = number + '';
    while (output.length < targetLength) {
      output = '0' + output;
    }
    return output;
  },
  numworks_maxtarsize:0x600000-0x10000,
  arrayBufferTransfer:function(oldBuffer, newByteLength) {
    const srcArray  = new Uint8Array(oldBuffer),
          destArray = new Uint8Array(newByteLength),
          copylen = Math.min(srcArray.buffer.byteLength, destArray.buffer.byteLength),
          floatArrayLength   = Math.trunc(copylen / 8),
          floatArraySource   = new Float64Array(srcArray.buffer,0,floatArrayLength),
          floarArrayDest     = new Float64Array(destArray.buffer,0,floatArrayLength);
    floarArrayDest.set(floatArraySource);
    let bytesCopied = floatArrayLength * 8;
    // slowpoke copy up to 7 bytes.
    while (bytesCopied < copylen ) {
      destArray[bytesCopied]=srcArray[bytesCopied];
      bytesCopied++;
    }
    return destArray.buffer;
  },
  tar_filesize:function(s){
    return 512*(1+Math.ceil(s/512));
  },
  _readString:function(buffer,str_offset, size) {
    let strView = new Uint8Array(buffer, str_offset, size);
    let i = 0;
    let rtnStr = "";
    while(i<size) {
      let ch=strView[i];
      if (ch<32) break;
      rtnStr += String.fromCharCode(ch);
      i++;
    }
    return rtnStr;
  },
  _readFileName:function(buffer,header_offset) {
    return UI._readString(buffer,header_offset, 100);
  },
  _readFileType:function(buffer,header_offset) {
    // offset: 156
    let typeView = new Uint8Array(buffer, header_offset+156, 1);
    let typeStr = String.fromCharCode(typeView[0]);
    if(typeStr == "0") {
      return "file";
    } else if(typeStr == "5") {
      return "directory";
    } else {
      return typeStr;
    }
  },
  _readFileSize:function(buffer,header_offset) {
    // offset: 124
    let szView = new Uint8Array(buffer, header_offset+124, 12);
    let szStr = "";
    for(let i = 0; i < 11; i++) {
      let tmp=szView[i];
      if (tmp<48 || tmp>57) return -1; // invalid file size
      szStr += String.fromCharCode(tmp);
    }
    return parseInt(szStr,8);
  },
  tar_clear:function(buffer){
    let view = new Uint8Array(buffer);
    for (let i=0;i<1024;++i)
      view[i]=0;
    $id('listtarfiles').innerHTML="";
    $id('perso').style.display='inline';    
  },
  tar_fileinfo:function(buffer,dohtml=1){
    let fileInfo = [];
    let offset = 0;
    let file_size = 0;       
    let file_name = "";
    let file_type = null;
    while (offset < buffer.byteLength - 512) {
      file_name = UI._readFileName(buffer,offset); // file name
      if (file_name.length == 0) 
        break;
      file_type = UI._readFileType(buffer,offset);
      file_size = UI._readFileSize(buffer,offset);
      if (file_size<0)
	break;
      //console.log(offset,file_name,file_size);
      fileInfo.push({
        "name": file_name,
        "type": file_type,
        "size": file_size,
        "header_offset": offset
      });

      offset += UI.tar_filesize(file_size);
    }
    //console.log('fileinfo',offset,dohtml);
    if (dohtml) UI.fileinfo2html(buffer,fileInfo,'listtarfiles');
    return fileInfo;
  },
  is_py:function(name,ext='.py'){
    return name.length>3 && name.substr(name.length-3,3)==ext;
  },
  fileinfo2html:function(buffer,finfo,fieldname){
    let l=finfo.length,s="<ul style='max-height:300px;overflow:auto'>";
    for (let i=0;i<l;++i){
      let cur=finfo[i];
      let S="<li>";
      let name=cur.name;
      S += name;
      S += " (";
      S += cur.size;
      S += "): ";
      let is_py=UI.is_py(name),is_xw=UI.is_py(name,'.xw');
      if (is_py || is_xw){
	let view='';
	let offset=cur.header_offset+512;
	let buf=new Uint8Array(buffer);
	for (let j=0;j<cur.size;++j)
	  view += String.fromCharCode(buf[offset+j]);
	S += '<a href="xcas'+(UI.langue==-1?'fr':'en')+'.html#exec&python=4&'+(is_py?'py':'xw')+'=0,0,'+encodeURIComponent(view)+'" target="_blank">Test</a>, ';
      }
      S += "<button onclick='UI.tar_removefile(UI.numworks_buffer,";
      S += '"'+cur.name+'"';
      S += i<11?',2':',1';
      S +=");' title='Cliquer ici pour effacer ce fichier'>X</button>";
      S += "<button onclick='UI.tar_savefile(";
      S += '"'+cur.name+'"';
      S +=");' title='Cliquer ici pour sauvegarder ce fichier'>&#x1f4be;</button>";      
      s+=S;
    }
    s += "</ul>";
    if (UI.langue==-1)
      s += UI.xcasfrmsg;
    else
      s += UI.xcasenmsg;
    //console.log(s);
    // put it the list of tarfiles
    if (fieldname){
      let field=$id(fieldname);
      console.log(field);
      field.innerHTML=s;
      field.style.display='inline';
      field=field.parentNode.style.display='inline';
    }
    // $id('perso').style.display='inline';
    return s;
  },
  storage_records:[],
  store_savefile:function(filename){
    let rec,j=0;
    for (j=0;j<UI.storage_records.length;++j){
      rec=UI.storage_records[j];
      // console.log(filename,rec.name);
      if (rec.name==filename) break;
    }
    if (j==UI.storage_records.length) return;
    //console.log(filename,rec);
    if (rec.type=="py")
      saveAs( new Blob([rec.code], {type: "text/plain;charset=utf-8"}), rec.name+"."+rec.type);
    else {
      saveAs(rec.data, rec.name+"."+rec.type);
    }
  },
  store2html:function(finfo,fieldname){
    let l=finfo.length,s="<ul style='max-height:300px;overflow:auto'>";
    for (let i=0;i<l;++i){
      let cur=finfo[i];
      let S="<li>";
      S += cur.name+"."+cur.type;
      S += " (";
      S += (cur.type=="py")?(cur.data.size-1):cur.data.size;
      S += "): ";
      if (cur.type=="py")
	S += '<a href="xcas'+(UI.langue==-1?'fr':'en')+'.html#exec&python=4&py=0,0,'+encodeURIComponent(cur.code)+'" target="_blank">Test</a>, ';
      //S += "<button onclick='UI.store_removefile(UI.numworks_buffer,";
      //S += '"'+cur.name+'"';
      //S +=");' title='Cliquer ici pour effacer ce fichier'>X</button>";
      S += "<button onclick='UI.store_savefile(";
      S += '"'+cur.name+'"';
      S +=");' title='Cliquer ici pour sauvegarder ce fichier'>&#x1f4be;</button>";      
      s+=S;
    }
    s += "</ul>";
    if (UI.langue==-1)
      s += UI.xcasfrmsg;
    else
      s += UI.xcasenmsg;
    //console.log(s);
    // put it the list of tarfiles
    if (fieldname){
      s += '<button onclick="$id(\''+fieldname+'\').style.display=\'none\'">X</button>';
      let field=$id(fieldname);
      field.innerHTML=s;
      field.style.display='inline';
      return s;
    }
    // $id('perso').style.display='inline';
    return s;
  },
  tar_writestring:function(buffer,str, offset, size) {
    let strView = new Uint8Array(buffer, offset, size);
    for(let i = 0; i < size; i++) {
      if (i < str.length) {
        strView[i] = str.charCodeAt(i);
      } else {
        strView[i] = 0;
      }
    }
  },
  tar_writechecksum:function(buffer,header_offset) {
    // offset: 148
    UI.tar_writestring(buffer,"        ", header_offset+148, 8); // first fill with spaces

    // add up header bytes
    let header = new Uint8Array(buffer, header_offset, 512);
    let chksum = 0;
    for (let i = 0; i < 512; i++) {
      chksum += header[i];
    }
    UI.tar_writestring(buffer,UI.leftpad(chksum.toString(8),6), header_offset+148, 8);
    UI.tar_writestring(buffer," ",header_offset+155,1); // add space inside chksum field
  },
  tar_fillheader:function(buffer,offset,exec=0){
    let uid = 501;
    let gid = 20;
    let mode = exec?"755":"644";
    let mtime = Date.now();
    let user = "user";
    let group = "group";

    UI.tar_writestring(buffer,UI.leftpad(mode,6)+" ", offset+100, 8);  
    UI.tar_writestring(buffer,UI.leftpad(uid.toString(8),6)+" ",offset+108,8);
    UI.tar_writestring(buffer,UI.leftpad(gid.toString(8),6)+" ",offset+116,8);
    UI.tar_writestring(buffer,UI.leftpad(Math.trunc(mtime/1000).toString(8),11)+" ",offset+136,12);

    //UI.tar_writestring(buffer,"ustar", offset+257,6); // magic string
    //UI.tar_writestring(buffer,"00", offset+263,2); // magic version
    UI.tar_writestring(buffer,"ustar  ", offset+257,8);
    
    UI.tar_writestring(buffer,user, offset+265,32); // user
    UI.tar_writestring(buffer,group, offset+297,32); //group
    UI.tar_writestring(buffer,"000000 ",offset+329,7); //devmajor
    UI.tar_writestring(buffer,"000000 ",offset+337,7); //devmajor
    UI.tar_writechecksum(buffer,offset);
  },
  tar_adddata:function(buffer,filename,data,exec=0,updatehtml=1){
    console.log('tar_adddata',exec,UI.tar_first_modif_offset);
    let finfo=UI.tar_fileinfo(buffer,0);
    let s=finfo.length,offset=0;
    if (s!=0){
      let last=finfo[s-1];
      offset=last.header_offset;
      offset += UI.tar_filesize(last.size);
    }
    let newsize=offset+1024+data.byteLength;
    newsize=10240*Math.ceil(newsize/10240);
    if (newsize>UI.numworks_maxtarsize) return 0;
    // console.log(buffer.byteLength,newsize);
    // resize buffer
    if (buffer.byteLength<newsize){
      buffer=UI.arrayBufferTransfer(buffer,newsize);
      if (updatehtml) UI.numworks_buffer=buffer;
    }
    console.log(buffer.byteLength,newsize);
    let view = new Uint8Array(buffer);    
    // add header
    // fill header with 0
    for (let i=0;i<1024;++i)
      view[offset+i]=0;
    UI.tar_writestring(buffer,filename,offset,100); // filename
    UI.tar_writestring(buffer,UI.leftpad(data.byteLength.toString(8),11)+" ",offset+124,12);  // filesize
    UI.tar_writestring(buffer,"0",offset+156,1); // file type
    UI.tar_fillheader(buffer,offset,exec);
    //let tmp=new Uint8Array(buffer,offset); console.log(tmp);
    //console.log(data.byteLength);
    // copy data 
    let srcview = new Uint8Array(data, 0, data.byteLength);    
    for (let i=0;i<data.byteLength;++i)
      view[offset+512+i]=srcview[i];
    let F=UI.tar_fileinfo(buffer,updatehtml);
    // console.log('tar_adddata',F);
    return buffer;
  },
  tar_removefile:function(buffer,filename,really=1){
    let finfo=UI.tar_fileinfo(buffer,0);
    let s=finfo.length;
    if (s==0) return 0;
    let i = finfo.findIndex(info => info.name == filename);
    if (i<0 || i>=s) return 0;
    if (really){
      msg = (UI.langue==-1?'Effacer ':'Erase ')+filename+'. ';
      msg += UI.langue==-1?'Etes-vous sur?':'Are you sure?';
      if (really==2)
	msg += (UI.langue==-1?'\nSi vous confirmez, votre firmware ne pourra plus etre certifie!':'If you confirm, your firmware can not be certified!');
      if (!confirm(msg))
	return 0;
    }
    let info = finfo[i];
    let view = new Uint8Array(buffer);
    // move info.header_offset+info.size+512 to info.header_offset
    let target=info.header_offset;
    if (target<UI.tar_first_modif_offset){
      UI.tar_first_modif_offset=target;
      console.log('tar_removefile',UI.tar_first_modif_offset);
    }
    let src=target+UI.tar_filesize(info.size);
    let infoend=finfo[s-1];
    let end=infoend.header_offset+UI.tar_filesize(infoend.size);
    for (src;src<end;++src,++target)
      view[target]=view[src];
    for (;target<end;++target) // clear space after new end
      view[target]=0;
    return UI.tar_fileinfo(buffer,1);
  },
  tar_addfile:function(file){
    let reader = new FileReader();
    reader.onerror = function (evt) { }
    reader.readAsArrayBuffer(file);
    reader.onload = function(evt){
      let fname=file.name,exec=1;
      if (fname.length>4 && fname.substr(fname.length-4,4)==".tar"){
	let reader = new FileReader();
	reader.onerror = function (evt) { }
	reader.readAsArrayBuffer(file);
	reader.onload = function(evt){
	  let buffer=evt.target.result;
	  let finfo=UI.tar_fileinfo(buffer,1),l=finfo.length;
	  for (let i=0;i<l;++i){ 
	    let cur=finfo[i];
	    let name=cur.name;
	    // check if name should be archived (not archived if it is in apps)
	    let size=cur.size;
	    let offset=cur.header_offset+512;
	    let data=new Uint8Array(buffer,offset,size);
	    if (UI.tar_adddata(UI.numworks_buffer,name,data,0,1)==0){
	      alert(UI.langue==-1?'Pas d\'espace en memoire flash':'No space left in flash');
	      return;
	    }
	  }
	}
	return;
      }
      if (fname.length>4 && fname.substr(fname.length-4,4)==".zip"){
	let zip=new JSZip();
	zip.loadAsync(evt.target.result).then(
	  function (zip) {
	    let lst=Object.entries(zip.files);
	    console.log(zip.files,lst);
	    for (let i=0;i<lst.length;++i){
	      let cur=lst[i];
	      let curname=cur[1].name;
	      console.log(curname);
	      zip.file(curname).async("arraybuffer").then(
		function(data){
		  exec=1;
		  for (let j=0;j<curname.length;++j){
		    if (curname[j]=='.'){
		      exec=0; break;
		    }
		  }
		  if (UI.tar_adddata(UI.numworks_buffer,curname,data,exec)==0){
		    alert(UI.langue==-1?'Pas d\'espace en memoire flash':'No space left in flash');
		    return;
		  }
		}
	      );
	    }
	  }
	);
	return;
      }
      for (let i=0;i<fname.length;++i){
	if (fname[i]=='.'){
	  exec=0; break;
	}
      }
      console.log(fname,exec);
      if (UI.tar_adddata(UI.numworks_buffer,file.name,evt.target.result,exec)==0)
	alert(UI.langue==-1?'Pas d\'espace en memoire flash':'No space left in flash');
    }
  },
  file_gettar:function(file){
    console.log(file);
    if (UI.numworks_buffer!=0){
      if (!confirm('Attention, la personnalisation courante va etre effacee! Continuer?'))
	return;
    }
    let reader = new FileReader();
    reader.onerror = function (evt) { }
    reader.readAsArrayBuffer(file);
    reader.onload = function(evt){
      UI.tar_first_modif_offset=0;
      UI.numworks_buffer=evt.target.result;
      console.log(UI.tar_fileinfo(UI.numworks_buffer,1));
    }
  },
  file_savetar:function(buffer,filename){
    let finfo=UI.tar_fileinfo(buffer,0);
    let s=finfo.length;
    if (s==0) return 0;
    let infoend=finfo[s-1];
    let end=infoend.header_offset+UI.tar_filesize(infoend.size);
    // add at least 1024 bytes, set them to 0
    let newsize=end+1024;  //newsize=10240*Math.ceil(newsize/10240);
    if (buffer.byteLength<newsize)
      buffer=UI.numworks_buffer=UI.arrayBufferTransfer(buffer,newsize);
    let view=new Uint8Array(buffer,0,newsize); 
    for (let i=end;i<newsize;++i)
      view[i]=0;
    let blob = new Blob([view], {type: "octet/stream"});
    saveAs(blob, filename);  
  },
  tar_savefile:function(filename){
    let buffer=UI.numworks_buffer;
    let finfo=UI.tar_fileinfo(buffer,0);
    let s=finfo.length;
    if (s==0) return 0;
    for (let i=0;i<s;++i){
      let cur=finfo[i];
      if (cur.name!=filename) continue;
      let view=new Uint8Array(cur.size);
      let offset=cur.header_offset+512;
      let buf=new Uint8Array(buffer);
      for (let j=0;j<cur.size;++j)
	view[j]=buf[offset+j];
      //console.log(cur,offset,view);
      let blob = new Blob([view], {type: "octet/stream"});
      saveAs(blob, filename);  
      return;
    }
  },
  numworks_gettar:function(calc){
    if (calc<=0){
      UI.numworks_gettar_(calc);
      return;
    }
    UI.calc=2;
    UI.nws_connect();
    window.setTimeout(UI.numworks_gettar_,100,calc);
  },
  numworks_gettar_:async function(calc){
    console.log(calc,UI.calculator);
    if (calc>0 && (UI.calculator==0 || !UI.calculator_connected)){
      alert(UI.langue==-1?UI.chkmsgfr:chkmsgen);
      if (UI.calculator) UI.calculator.stopAutoConnect();
      return;
    }
    if (UI.numworks_buffer!=0){
      if (!confirm('Attention, la personnalisation courante va etre effacee! Continuer?'))
	return;
    }
    UI.tar_first_modif_offset=0;
    if (calc==0)
      UI.numworks_buffer=await UI.loadfile('apps.tar');
    else if (calc==-1)
      UI.numworks_buffer=await UI.loadfile('appsalpha.tar');
    else
      UI.numworks_buffer=await UI.calculator.get_apps();
    let finfo=UI.tar_fileinfo(UI.numworks_buffer,1),s=finfo.length;
    if (calc>0 && s>0)
      UI.tar_first_modif_offset=finfo[s-1].header_offset+UI.tar_filesize(finfo[s-1].size);
    // window.setTimeout(UI.tar_fileinfo,1000,UI.numworks_buffer);
    //UI.print("Archive flash:\n"+UI.numworks_tarinfo(UI.numworks_buffer));
  },
  numworks_getarchive:function(filename,appendapps=1){
    // filename=string filename to archive or 0 or 1 for calculator upgrade
    UI.calc=2;
    UI.nws_connect();
    window.setTimeout(UI.numworks_getarchive_,100,filename,appendapps);
  },
  numworks_getarchive_:async function(filename,appendapps=1){
    let sendcalc= !(typeof filename=="string");
    if (sendcalc) appendapps=1;
    console.log(sendcalc,UI.calculator);
    if (UI.calculator==0 || !UI.calculator_connected){
      alert(UI.langue==-1?UI.chkmsgfr:chkmsgen);
      if (UI.calculator) UI.calculator.stopAutoConnect();
      return;
    }
    let buffer=await UI.calculator.get_apps();
    let appsbuffer;
    if (filename==1)
      appsbuffer=await UI.loadfile('appsalpha.tar');
    else
      appsbuffer=await UI.loadfile('apps.tar');
    let appsfinfo=UI.tar_fileinfo(appsbuffer);
    console.log('appsfinfo',appsfinfo);
    let resbuffer;
    if (appendapps)
      resbuffer=appsbuffer;
    else
      resbuffer=new Uint8Array().buffer;
    let finfo=UI.tar_fileinfo(buffer,0),l=finfo.length,L=appsfinfo.length;
    for (let i=0;i<l;++i){ 
      let cur=finfo[i];
      let name=cur.name,j;
      // check if name should be archived (not archived if it is in apps)
      for (j=0;j<L;++j){
	if (name==appsfinfo[j].name)
	  break;
      }
      console.log(name,j,L);
      if (j<L)
	continue;
      let size=cur.size;
      let offset=cur.header_offset+512;
      let data=new Uint8Array(buffer,offset,size);
      //console.log(data);
      resbuffer=UI.tar_adddata(resbuffer,name,data,0,0);
      // console.log('getarchive',i,UI.tar_fileinfo(resbuffer,0));
    }
    finfo=UI.tar_fileinfo(resbuffer,0);
    console.log(filename,finfo);
    //return;
    if (sendcalc){
      UI.calculator.device.startAddress = 0x90200000;
      UI.calculator.device.logProgress(0,'apps');
      let res=await UI.calculator.device.do_download(UI.calculator.transferSize, resbuffer, false);
      alert(UI.langue==-1?'Installation terminee. Appuyer sur RESET a l\'arriere de la calculatrice':'Install success. Press RESET on the calculator back.');
      return res;
    }
    else
      UI.file_savetar(resbuffer,filename);
  },
  numworks_sendtar:function(){
    UI.calc=2;
    UI.nws_connect();
    window.setTimeout(UI.numworks_sendtar_,100);
  },
  numworks_sendtar_:async function(){
    if (!UI.numworks_chk_buffer())
      return;
    if (!UI.calculator || !UI.calculator_connected){
      alert(UI.langue==-1?UI.chkmsgfr:chkmsgen);
      if (UI.calculator) UI.calculator.stopAutoConnect();
      return;
    }
    $id('perso').style.display='none';    
    if (UI.tar_first_modif_offset){
      let finfo=UI.tar_fileinfo(UI.numworks_buffer,0);
      let s=finfo.length;
      if (s==0) return;
      let start=65536*Math.floor(UI.tar_first_modif_offset/65536);
      let size=finfo[s-1].header_offset+UI.tar_filesize(finfo[s-1].size)-start;
      let buf=new Uint8Array(size);
      let src=new Uint8Array(UI.numworks_buffer);
      for (let i=0;i<size;++i)
	buf[i]=src[start+i];
      let buffer=buf.buffer;
      console.log('sendtar',UI.tar_first_modif_offset,start,size,buffer.byteLength);
      UI.calculator.device.startAddress = 0x90200000+start;
      UI.calculator.device.logProgress(0,'apps');
      res=await UI.calculator.device.do_download(UI.calculator.transferSize, buffer, false);      
    }
    else {
      UI.calculator.device.startAddress = 0x90200000;
      UI.calculator.device.logProgress(0,'apps');
      res=await UI.calculator.device.do_download(UI.calculator.transferSize, UI.numworks_buffer, false);
    }
    $id('perso').style.display='inline';    
    UI.print('apps OK, press RESET on calculator back');
  },
  numworks_certify:function(sigfile,rwcheck=false){
    $id('persoinit').style.display='none';    
    UI.calc=2;
    UI.nws_connect();
    window.setTimeout(UI.numworks_certify_,100,sigfile,rwcheck);
  },
  numworks_certify_:async function(sigfile,rwcheck){
    if (UI.calculator==0 || !UI.calculator_connected){
      alert(UI.langue==-1?UI.chkmsgfr:chkmsgen);
      if (UI.calculator) UI.calculator.stopAutoConnect();
      return -1;
    }
    if (rwcheck)
      if (!confirm(UI.langue==-1?"Ce test necessite l'accord du proprietaire de la calculatrice, il dure environ 1 minute. Poursuivre?":'This test requires the calculator owner authorization. It will take about 1 minute. Perform?'))
	return -1;
    //else alert(UI.langue==-1?'Le test va prendre environ 20 secondes':'Test will take about 20 seconds');
    UI.print('========');
    let res;
    let external=await UI.calculator.get_external_flash(1);
    UI.adjust_exam_sector(external);
    res=await UI.sig_check(sigfile,external,'delta.external.bin');
    if (!res){
      alert(UI.langue==-1?'Flash externe slot 1 non certifiee. Cliquer sur installer KhiCAS pour installer un firmware certifie.':'External flash slot 1 not certified');
      return 2;
    }
    UI.print('External flash slot 1 OK');
    external=await UI.calculator.get_external_flash(2);
    UI.adjust_exam_sector(external);
    res=await UI.sig_check(sigfile,external,'khi.B.bin');
    if (!res){
      alert(UI.langue==-1?'Flash externe slot 2 non certifiee. Cliquer sur installer KhiCAS pour installer un firmware certifie.':'External flash slot 2 not certified');
      return 2;
    }
    UI.print('External flash slot 2 OK');
    let apps=await UI.calculator.get_apps();
    res=await UI.sig_check(sigfile,apps,'apps.tar');
    if (!res){
      alert(UI.langue==-1?'Applications non certifiees. Cliquer sur installer KhiCAS pour installer un firmware certifie.':'Applications not certified');
      return 3;
    }
    UI.print('Apps OK');
    if (rwcheck){
      if (0){
	UI.print('R/W check 0x91300000');
	res=await UI.calculator.rw_check(0x90130000,0xd0000);
	if (!res){
	  alert(UI.langue==-1?'Echec du test lecture/ecriture':'Read/Write test failure');
	  return 4;
	}
      }
      UI.print('R/W check 0x90740000');
      res=await UI.calculator.rw_check(0x90740000,0xa0000);
      if (!res){
	alert(UI.langue==-1?'Echec du test lecture/ecriture':'Read/Write test failure');
	return 4;
      }
      UI.print('R/W OK');
    }
    let pinfo = await UI.calculator.getPlatformInfo();
    if (pinfo["bootloader"]==0){
      let internal=await UI.calculator.get_internal_flash();
      //console.log(sigfile);
      res=await UI.sig_check(sigfile,internal,'delta.internal.bin');
      if (!res){
	alert(UI.langue==-1?'Flash interne non certifiee. Cliquer sur installer KhiCAS pour installer un firmware certifie.':'Internal flash not certified');
	return 1;
      }
      UI.print('Internal flash OK');
      alert(UI.langue==-1?'Firmware certifie (interne, externe et KhiCAS)':'Certified firmware (internal, external and KhiCAS)');
    }
    else
      alert(UI.langue==-1?'Firmware certifie (externe et KhiCAS). Pour certifier la flash interne, faire Reset et lancer KhiCAS.':'Firmware certified (external and KhiCAS). Reset the calculator and run KhiCAS to certify internal flash.');      
    return 0;
  },
  htmlbuffer:'',
  htmlcheck:true,
  print:function(text){
    let element = $id('output');
    //console.log(text.charCodeAt(0));
    if (text.length == 1 && text.charCodeAt(0) == 12) {
      element.innerHTML = '';
      return;
    }
    if (text.length >= 1 && text.charCodeAt(0) == 2) {
      console.log('STX');
      UI.htmlcheck = false;
      UI.htmlbuffer = '';
      return;
    }
    if (text.length >= 1 && text.charCodeAt(0) == 3) {
      console.log('ETX');
      UI.htmlcheck = true;
      element.style.display = 'inherit';
      element.innerHTML += UI.htmlbuffer;
      UI.htmlbuffer = '';
      element.scrollTop = 99999;
      return;
    }
    if (UI.htmlcheck) {
      // These replacements are necessary if you render to raw HTML
      text = '' + text;
      console.log(text);
      text = text.replace(/&/g, "&amp;");
      text = text.replace(/</g, "&lt;");
      text = text.replace(/>/g, "&gt;");
      text = text.replace(/\n/g, '<br>');
      text += '<br>';
      element.style.display = 'inherit';
      element.innerHTML += text; // element.value += text + "\n";
      element.scrollTop = 99999; // focus on bottom
    } else UI.htmlbuffer += text;
    element.scrollIntoView();
  },
};

