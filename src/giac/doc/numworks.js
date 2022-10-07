// (c) Maxime Friess, 2021, licence MIT, https://github.com/M4xi1m3/numworks.js/blob/master/README.md
var Numworks =
/******/ (function(modules) { // webpackBootstrap
/******/ 	// The module cache
/******/ 	var installedModules = {};
/******/
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/
/******/ 		// Check if module is in cache
/******/ 		if(installedModules[moduleId]) {
/******/ 			return installedModules[moduleId].exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = installedModules[moduleId] = {
/******/ 			i: moduleId,
/******/ 			l: false,
/******/ 			exports: {}
/******/ 		};
/******/
/******/ 		// Execute the module function
/******/ 		modules[moduleId].call(module.exports, module, module.exports, __webpack_require__);
/******/
/******/ 		// Flag the module as loaded
/******/ 		module.l = true;
/******/
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/
/******/
/******/ 	// expose the modules object (__webpack_modules__)
/******/ 	__webpack_require__.m = modules;
/******/
/******/ 	// expose the module cache
/******/ 	__webpack_require__.c = installedModules;
/******/
/******/ 	// define getter function for harmony exports
/******/ 	__webpack_require__.d = function(exports, name, getter) {
/******/ 		if(!__webpack_require__.o(exports, name)) {
/******/ 			Object.defineProperty(exports, name, { enumerable: true, get: getter });
/******/ 		}
/******/ 	};
/******/
/******/ 	// define __esModule on exports
/******/ 	__webpack_require__.r = function(exports) {
/******/ 		if(typeof Symbol !== 'undefined' && Symbol.toStringTag) {
/******/ 			Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
/******/ 		}
/******/ 		Object.defineProperty(exports, '__esModule', { value: true });
/******/ 	};
/******/
/******/ 	// create a fake namespace object
/******/ 	// mode & 1: value is a module id, require it
/******/ 	// mode & 2: merge all properties of value into the ns
/******/ 	// mode & 4: return value when already ns object
/******/ 	// mode & 8|1: behave like require
/******/ 	__webpack_require__.t = function(value, mode) {
/******/ 		if(mode & 1) value = __webpack_require__(value);
/******/ 		if(mode & 8) return value;
/******/ 		if((mode & 4) && typeof value === 'object' && value && value.__esModule) return value;
/******/ 		var ns = Object.create(null);
/******/ 		__webpack_require__.r(ns);
/******/ 		Object.defineProperty(ns, 'default', { enumerable: true, value: value });
/******/ 		if(mode & 2 && typeof value != 'string') for(var key in value) __webpack_require__.d(ns, key, function(key) { return value[key]; }.bind(null, key));
/******/ 		return ns;
/******/ 	};
/******/
/******/ 	// getDefaultExport function for compatibility with non-harmony modules
/******/ 	__webpack_require__.n = function(module) {
/******/ 		var getter = module && module.__esModule ?
/******/ 			function getDefault() { return module['default']; } :
/******/ 			function getModuleExports() { return module; };
/******/ 		__webpack_require__.d(getter, 'a', getter);
/******/ 		return getter;
/******/ 	};
/******/
/******/ 	// Object.prototype.hasOwnProperty.call
/******/ 	__webpack_require__.o = function(object, property) { return Object.prototype.hasOwnProperty.call(object, property); };
/******/
/******/ 	// __webpack_public_path__
/******/ 	__webpack_require__.p = "";
/******/
/******/
/******/ 	// Load entry module and return exports
/******/ 	return __webpack_require__(__webpack_require__.s = 3);
/******/ })
/************************************************************************/
/******/ ([
/* 0 */
/***/ (function(module, exports, __webpack_require__) {


const DFU = __webpack_require__(1);
const DFUse = __webpack_require__(5);

module.exports = {
    DFU: DFU,
    DFUse: DFUse
};


/***/ }),
/* 1 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


/*
 * Static DFU class.
 */
class DFU {
    static get DETACH() { return 0x00; }
    static get DNLOAD() { return 0x01; }
    static get UPLOAD() { return 0x02; }
    static get GETSTATUS() { return 0x03; }
    static get CLRSTATUS() { return 0x04; }
    static get GETSTATE() { return 0x05; }
    static get ABORT() { return 0x06; }
    
    static get appIDLE() { return 0; }
    static get appDETACH() { return 1; }
    static get dfuIDLE() { return 2; }
    static get dfuDNLOAD_SYNC() { return 3; }
    static get dfuDNBUSY() { return 4; }
    static get dfuDNLOAD_IDLE() { return 5; }
    static get dfuMANIFEST_SYNC() { return 6; }
    static get dfuMANIFEST() { return 7; }
    static get dfuMANIFEST_WAIT_RESET() { return 8; }
    static get dfuUPLOAD_IDLE() { return 9; }
    static get dfuERROR() { return 10; }
    
    static get STATUS_OK() { return 0x0; }
    
    // Device = null;
    
    static findDeviceDfuInterfaces(device) {
        let interfaces = [];
        for (let conf of device.configurations) {
            for (let intf of conf.interfaces) {
                for (let alt of intf.alternates) {
                    if (alt.interfaceClass === 0xFE &&
                        alt.interfaceSubclass === 0x01 &&
                        (alt.interfaceProtocol === 0x01 || alt.interfaceProtocol === 0x02)) {
                        let settings = {
                            "configuration": conf,
                            "interface": intf,
                            "alternate": alt,
                            "name": alt.interfaceName
                        };
                        interfaces.push(settings);
                    }
                }
            }
        }

        return interfaces;
    }
    
    static findAllDfuInterfaces() {
        return navigator.usb.getDevices().then(
            devices => {
                let matches = [];
                for (let device of devices) {
                    let interfaces = DFU.findDeviceDfuInterfaces(device);
                    for (let interface_ of interfaces) {
                        matches.push(new DFU.Device(device, interface_));
                    }
                }
                return matches;
            }
        )
    }
    
    static parseDeviceDescriptor(data) {
        return {
            bLength:            data.getUint8(0),
            bDescriptorType:    data.getUint8(1),
            bcdUSB:             data.getUint16(2, true),
            bDeviceClass:       data.getUint8(4),
            bDeviceSubClass:    data.getUint8(5),
            bDeviceProtocol:    data.getUint8(6),
            bMaxPacketSize:     data.getUint8(7),
            idVendor:           data.getUint16(8, true),
            idProduct:          data.getUint16(10, true),
            bcdDevice:          data.getUint16(12, true),
            iManufacturer:      data.getUint8(14),
            iProduct:           data.getUint8(15),
            iSerialNumber:      data.getUint8(16),
            bNumConfigurations: data.getUint8(17),
        };
    }
    
    static parseConfigurationDescriptor(data) {
        let descriptorData = new DataView(data.buffer.slice(9));
        let descriptors = DFU.parseSubDescriptors(descriptorData);
        return {
            bLength:            data.getUint8(0),
            bDescriptorType:    data.getUint8(1),
            wTotalLength:       data.getUint16(2, true),
            bNumInterfaces:     data.getUint8(4),
            bConfigurationValue:data.getUint8(5),
            iConfiguration:     data.getUint8(6),
            bmAttributes:       data.getUint8(7),
            bMaxPower:          data.getUint8(8),
            descriptors:        descriptors
        };
    }

    static parseInterfaceDescriptor(data) {
        return {
            bLength:            data.getUint8(0),
            bDescriptorType:    data.getUint8(1),
            bInterfaceNumber:   data.getUint8(2),
            bAlternateSetting:  data.getUint8(3),
            bNumEndpoints:      data.getUint8(4),
            bInterfaceClass:    data.getUint8(5),
            bInterfaceSubClass: data.getUint8(6),
            bInterfaceProtocol: data.getUint8(7),
            iInterface:         data.getUint8(8),
            descriptors:        []
        };
    }

    static parseFunctionalDescriptor(data) {
        return {
            bLength:           data.getUint8(0),
            bDescriptorType:   data.getUint8(1),
            bmAttributes:      data.getUint8(2),
            wDetachTimeOut:    data.getUint16(3, true),
            wTransferSize:     data.getUint16(5, true),
            bcdDFUVersion:     data.getUint16(7, true)
        };
    }

    static parseSubDescriptors(descriptorData) {
        const DT_INTERFACE = 4;
        // const DT_ENDPOINT = 5;
        const DT_DFU_FUNCTIONAL = 0x21;
        const USB_CLASS_APP_SPECIFIC = 0xFE;
        const USB_SUBCLASS_DFU = 0x01;
        let remainingData = descriptorData;
        let descriptors = [];
        let currIntf;
        let inDfuIntf = false;
        while (remainingData.byteLength > 2) {
            let bLength = remainingData.getUint8(0);
            let bDescriptorType = remainingData.getUint8(1);
            let descData = new DataView(remainingData.buffer.slice(0, bLength));
            if (bDescriptorType === DT_INTERFACE) {
                currIntf = DFU.parseInterfaceDescriptor(descData);
                if (currIntf.bInterfaceClass === USB_CLASS_APP_SPECIFIC &&
                    currIntf.bInterfaceSubClass === USB_SUBCLASS_DFU) {
                    inDfuIntf = true;
                } else {
                    inDfuIntf = false;
                }
                descriptors.push(currIntf);
            } else if (inDfuIntf && bDescriptorType === DT_DFU_FUNCTIONAL) {
                let funcDesc = DFU.parseFunctionalDescriptor(descData)
                descriptors.push(funcDesc);
                currIntf.descriptors.push(funcDesc);
            } else {
                let desc = {
                    bLength: bLength,
                    bDescriptorType: bDescriptorType,
                    data: descData
                };
                descriptors.push(desc);
                if (currIntf) {
                    currIntf.descriptors.push(desc);
                }
            }
            remainingData = new DataView(remainingData.buffer.slice(bLength));
        }

        return descriptors;
    }
}

/**
 * Represents a DFU-enabled connected device.
 */
DFU.Device = class {
    constructor(device, settings) {
        this.device_ = device;
        this.settings = settings;
        this.intfNumber = settings["interface"].interfaceNumber;
        this.dnload = this.download;
        this.clrStatus = this.clearStatus;
    }
    
    logDebug(msg) {
        console.debug(msg);
    }
    
    logInfo(msg) {
        console.info(msg);
    }
    
    logWarning(msg) {
        console.warn(msg);
    }
    
    logError(msg) {
        console.error(msg);
    }
    
  logProgress(done, total) {
    //console.log('logprogress',UI.nws_progress,done,total);
    if (typeof UI!=='undefined' && typeof UI.nws_progress!=='undefined'){
      UI.nws_progress.style.display='inline';
      //let mainview=document.getElementById('report');
      //if (typeof mainview!=='undefined') mainview.scrollIntoView();
      if (typeof total !== 'undefined'){
	if (typeof total=='string')
	  UI.nws_progresslegend.innerHTML=total;
	else {
	  UI.nws_progress.max=total;
	  if (done>=total){
	    UI.nws_progress.style.display='none';
	    UI.nws_progresslegend.innerHTML='';
	    //mainview=document.getElementById('giac');
	    //if (typeof mainview!=='undefined') mainview.scrollIntoView();	    
	  }
	}
      }
      UI.nws_progress.value=done;
      return;
    }
    if (typeof total === 'undefined') {
      this.logDebug(done)
    } else {
      this.logDebug(done + '/' + total);
    }
  }
    
    async open() {
        await this.device_.open();
        const confValue = this.settings.configuration.configurationValue;
        if (this.device_.configuration === null ||
            this.device_.configuration.configurationValue !== confValue) {
            await this.device_.selectConfiguration(confValue);
        }

        const intfNumber = this.settings["interface"].interfaceNumber;
        if (!this.device_.configuration.interfaces[intfNumber].claimed) {
            await this.device_.claimInterface(intfNumber);
        }

        const altSetting = this.settings.alternate.alternateSetting;
        let intf = this.device_.configuration.interfaces[intfNumber];
        if (intf.alternate === null ||
            intf.alternate.alternateSetting !== altSetting) {
            await this.device_.selectAlternateInterface(intfNumber, altSetting);
        }
    }
    
    async close() {
        try {
            await this.device_.close();
        } catch (error) {
            console.log(error);
        }
    }
    
    readDeviceDescriptor() {
        const GET_DESCRIPTOR = 0x06;
        const DT_DEVICE = 0x01;
        const wValue = (DT_DEVICE << 8);

        return this.device_.controlTransferIn({
            "requestType": "standard",
            "recipient": "device",
            "request": GET_DESCRIPTOR,
            "value": wValue,
            "index": 0
        }, 18).then(
            result => {
                if (result.status === "ok") {
                     return Promise.resolve(result.data);
                } else {
                    return Promise.reject(result.status);
                }
            }
        );
    }
    
    async readStringDescriptor(index, langID) {
        if (typeof langID === 'undefined') {
            langID = 0;
        }

        const GET_DESCRIPTOR = 0x06;
        const DT_STRING = 0x03;
        const wValue = (DT_STRING << 8) | index;

        const request_setup = {
            "requestType": "standard",
            "recipient": "device",
            "request": GET_DESCRIPTOR,
            "value": wValue,
            "index": langID
        }

        // Read enough for bLength
        var result = await this.device_.controlTransferIn(request_setup, 1);

        if (result.status === "ok") {
            // Retrieve the full descriptor
            const bLength = result.data.getUint8(0);
            result = await this.device_.controlTransferIn(request_setup, bLength);
            if (result.status === "ok") {
                const len = (bLength-2) / 2;
                let u16_words = [];
                for (let i=0; i < len; i++) {
                    u16_words.push(result.data.getUint16(2+i*2, true));
                }
                if (langID === 0) {
                    // Return the langID array
                    return u16_words;
                } else {
                    // Decode from UCS-2 into a string
                    return String.fromCharCode.apply(String, u16_words);
                }
            }
        }
        
        throw new Error(`Failed to read string descriptor ${index}: ${result.status}`);
    }
    
    async readInterfaceNames() {
        const DT_INTERFACE = 4;

        let configs = {};
        let allStringIndices = new Set();
        for (let configIndex=0; configIndex < this.device_.configurations.length; configIndex++) {
            const rawConfig = await this.readConfigurationDescriptor(configIndex);
            let configDesc = DFU.parseConfigurationDescriptor(rawConfig);
            let configValue = configDesc.bConfigurationValue;
            configs[configValue] = {};

            // Retrieve string indices for interface names
            for (let desc of configDesc.descriptors) {
                if (desc.bDescriptorType === DT_INTERFACE) {
                    if (!(desc.bInterfaceNumber in configs[configValue])) {
                        configs[configValue][desc.bInterfaceNumber] = {};
                    }
                    configs[configValue][desc.bInterfaceNumber][desc.bAlternateSetting] = desc.iInterface;
                    if (desc.iInterface > 0) {
                        allStringIndices.add(desc.iInterface);
                    }
                }
            }
        }

        let strings = {};
        // Retrieve interface name strings
        for (let index of allStringIndices) {
            try {
                strings[index] = await this.readStringDescriptor(index, 0x0409);
            } catch (error) {
                console.log(error);
                strings[index] = null;
            }
        }

        for (let configValue in configs) {
            for (let intfNumber in configs[configValue]) {
                for (let alt in configs[configValue][intfNumber]) {
                    const iIndex = configs[configValue][intfNumber][alt];
                    configs[configValue][intfNumber][alt] = strings[iIndex];
                }
            }
        }

        return configs;
    }
    
    readConfigurationDescriptor(index) {
        const GET_DESCRIPTOR = 0x06;
        const DT_CONFIGURATION = 0x02;
        const wValue = ((DT_CONFIGURATION << 8) | index);

        return this.device_.controlTransferIn({
            "requestType": "standard",
            "recipient": "device",
            "request": GET_DESCRIPTOR,
            "value": wValue,
            "index": 0
        }, 4).then(
            result => {
                if (result.status === "ok") {
                    // Read out length of the configuration descriptor
                    let wLength = result.data.getUint16(2, true);
                    return this.device_.controlTransferIn({
                        "requestType": "standard",
                        "recipient": "device",
                        "request": GET_DESCRIPTOR,
                        "value": wValue,
                        "index": 0
                    }, wLength);
                } else {
                    return Promise.reject(result.status);
                }
            }
        ).then(
            result => {
                if (result.status === "ok") {
                    return Promise.resolve(result.data);
                } else {
                    return Promise.reject(result.status);
                }
            }
        );
    }
    
    requestOut(bRequest, data, wValue=0) {
        return this.device_.controlTransferOut({
            "requestType": "class",
            "recipient": "interface",
            "request": bRequest,
            "value": wValue,
            "index": this.intfNumber
        }, data).then(
            result => {
                if (result.status === "ok") {
                    return Promise.resolve(result.bytesWritten);
                } else {
                    return Promise.reject(result.status);
                }
            },
            error => {
                return Promise.reject("ControlTransferOut failed: " + error);
            }
        );
    }
    
    requestIn(bRequest, wLength, wValue=0) {
        return this.device_.controlTransferIn({
            "requestType": "class",
            "recipient": "interface",
            "request": bRequest,
            "value": wValue,
            "index": this.intfNumber
        }, wLength).then(
            result => {
                if (result.status === "ok") {
                    return Promise.resolve(result.data);
                } else {
                    return Promise.reject(result.status);
                }
            },
            error => {
                console.error(error);
                return Promise.reject("ControlTransferIn failed: " + error);
            }
        );
    }
    
    detach() {
        return this.requestOut(DFU.DETACH, undefined, 1000);
    }
    
    async waitDisconnected(timeout) {
        let device = this;
        let usbDevice = this.device_;
        return new Promise(function(resolve, reject) {
            let timeoutID;
            if (timeout > 0) {
                /*
                function onTimeout() {
                    navigator.usb.removeEventListener("disconnect", onDisconnect);
                    if (device.disconnected !== true) {
                        reject("Disconnect timeout expired");
                    }
                }
                */
                timeoutID = setTimeout(reject, timeout);
            }

            function onDisconnect(event) {
                if (event.device === usbDevice) {
                    if (timeout > 0) {
                        clearTimeout(timeoutID);
                    }
                    device.disconnected = true;
                    navigator.usb.removeEventListener("disconnect", onDisconnect);
                    event.stopPropagation();
                    resolve(device);
                }
            }

            navigator.usb.addEventListener("disconnect", onDisconnect);
        });
    }
    
    download(data, blockNum) {
        return this.requestOut(DFU.DNLOAD, data, blockNum);
    }
    
    upload(length, blockNum) {
        return this.requestIn(DFU.UPLOAD, length, blockNum)
    }
    
    clearStatus() {
        return this.requestOut(DFU.CLRSTATUS);
    }
    
    getStatus() {
        return this.requestIn(DFU.GETSTATUS, 6).then(
            data =>
                Promise.resolve({
                    "status": data.getUint8(0),
                    "pollTimeout": data.getUint32(1, true) & 0xFFFFFF,
                    "state": data.getUint8(4)
                }),
            error =>
                Promise.reject("DFU GETSTATUS failed: " + error)
        );
    }
    
    getState() {
        return this.requestIn(DFU.GETSTATE, 1).then(
            data => Promise.resolve(data.getUint8(0)),
            error => Promise.reject("DFU GETSTATE failed: " + error)
        );
    }
    
    abort() {
        return this.requestOut(DFU.ABORT);
    }
    
    async abortToIdle() {
        await this.abort();
        let state = await this.getState();
        if (state === DFU.dfuERROR) {
            await this.clearStatus();
            state = await this.getState();
        }
        if (state !== DFU.dfuIDLE) {
            throw new Error("Failed to return to idle state after abort: state " + state.state);
        }
    }
    
    async do_upload(xfer_size, max_size=Infinity, first_block=0) {
        let transaction = first_block;
        let blocks = [];
        let bytes_read = 0;

        this.logInfo("Copying data from DFU device to browser");
        // Initialize progress to 0
        this.logProgress(0);

        let result;
        let bytes_to_read;
        do {
            bytes_to_read = Math.min(xfer_size, max_size - bytes_read);
            result = await this.upload(bytes_to_read, transaction++);
            // this.logDebug("Read " + result.byteLength + " bytes");
            if (result.byteLength > 0) {
                blocks.push(result);
                bytes_read += result.byteLength;
            }
            if (Number.isFinite(max_size)) {
                this.logProgress(bytes_read, max_size);
            } else {
                this.logProgress(bytes_read);
            }
        } while ((bytes_read < max_size) && (result.byteLength === bytes_to_read));

        if (bytes_read === max_size) {
            await this.abortToIdle();
        }

        this.logInfo(`Read ${bytes_read} bytes`);

        return new Blob(blocks, { type: "application/octet-stream" });
    }
    
    async poll_until(state_predicate) {
        let dfu_status = await this.getStatus();

        let device = this;
        function async_sleep(duration_ms) {
            return new Promise(function(resolve, reject) {
                device.logDebug("Sleeping for " + duration_ms + "ms");
                setTimeout(resolve, duration_ms);
            });
        }
        
        while (!state_predicate(dfu_status.state) && dfu_status.state !== DFU.dfuERROR) {
            await async_sleep(dfu_status.pollTimeout);
            dfu_status = await this.getStatus();
        }

        return dfu_status;
    }
    
    poll_until_idle(idle_state) {
        return this.poll_until(state => (state === idle_state));
    }
    
    async do_download(xfer_size, data, manifestationTolerant) {
        let bytes_sent = 0;
        let expected_size = data.byteLength;
        let transaction = 0;

        this.logInfo("Copying data from browser to DFU device");

        // Initialize progress to 0
        this.logProgress(bytes_sent, expected_size);

        while (bytes_sent < expected_size) {
            const bytes_left = expected_size - bytes_sent;
            const chunk_size = Math.min(bytes_left, xfer_size);

            let bytes_written = 0;
            let dfu_status;
            try {
                bytes_written = await this.download(data.slice(bytes_sent, bytes_sent+chunk_size), transaction++);
                this.logDebug("Sent " + bytes_written + " bytes");
                dfu_status = await this.poll_until_idle(DFU.dfuDNLOAD_IDLE);
            } catch (error) {
                throw new Error("Error during DFU download: " + error);
            }

            if (dfu_status.status !== DFU.STATUS_OK) {
                throw new Error(`DFU DOWNLOAD failed state=${dfu_status.state}, status=${dfu_status.status}`);
            }

            this.logDebug("Wrote " + bytes_written + " bytes");
            bytes_sent += bytes_written;

            this.logProgress(bytes_sent, expected_size);
        }

        this.logDebug("Sending empty block");
        try {
            await this.download(new ArrayBuffer([]), transaction++);
        } catch (error) {
            throw new Error("Error during final DFU download: " + error);
        }

        this.logInfo("Wrote " + bytes_sent + " bytes");
        this.logInfo("Manifesting new firmware");

        if (manifestationTolerant) {
            // Transition to MANIFEST_SYNC state
            let dfu_status;
            try {
                // Wait until it returns to idle.
                // If it's not really manifestation tolerant, it might transition to MANIFEST_WAIT_RESET
                dfu_status = await this.poll_until(state => (state === DFU.dfuIDLE || state === DFU.dfuMANIFEST_WAIT_RESET));
                if (dfu_status.state === DFU.dfuMANIFEST_WAIT_RESET) {
                    this.logDebug("Device transitioned to MANIFEST_WAIT_RESET even though it is manifestation tolerant");
                }
                if (dfu_status.status !== DFU.STATUS_OK) {
                    throw new Error(`DFU MANIFEST failed state=${dfu_status.state}, status=${dfu_status.status}`);
                }
            } catch (error) {
                if (error.endsWith("ControlTransferIn failed: NotFoundError: Device unavailable.") ||
                    error.endsWith("ControlTransferIn failed: NotFoundError: The device was disconnected.")) {
                    this.logWarning("Unable to poll final manifestation status");
                } else {
                    throw new Error("Error during DFU manifest: " + error);
                }
            }
        } else {
            // Try polling once to initiate manifestation
            try {
                let final_status = await this.getStatus();
                this.logDebug(`Final DFU status: state=${final_status.state}, status=${final_status.status}`);
            } catch (error) {
                this.logDebug("Manifest GET_STATUS poll error: " + error);
            }
        }

        // Reset to exit MANIFEST_WAIT_RESET
        try {
            await this.device_.reset();
        } catch (error) {
            if (error === "NetworkError: Unable to reset the device." ||
                error === "NotFoundError: Device unavailable." ||
                error === "NotFoundError: The device was disconnected.") {
                this.logDebug("Ignored reset error");
            } else {
                throw new Error("Error during reset for manifestation: " + error);
            }
        }

        return;
    };
};

module.exports = DFU;



/***/ }),
/* 2 */
/***/ (function(module, exports) {


/**
 * Class to parse and reconstruct the numworks' internal storage.
 * Only parses python scripts for now, ditches the rest.
 * @TODO parse other things.
 *
 * @author Maxime "M4x1m3" FRIESS
 * @license MIT
 */
class Storage {
    constructor() {
        this.magik = null;
        this.records = null;
    }
    
    async __encodePyRecord(record) {
        var content = new TextEncoder("utf-8").encode(record.code);
        
        record.data = new Blob([
            concatTypedArrays(
                new Uint8Array([record.autoImport ? 1 : 0]),
                concatTypedArrays(
                    content,
                    new Uint8Array([0])
                )
            )
        ]);
        
        delete record.autoImport;
        delete record.code;
        
        return record;
    }
    
    __getRecordEncoders() {
        return {
            py: this.__encodePyRecord.bind(this)
        };
    }
    
    async __assembleStorage(records, maxSize) {
        const encoder = new TextEncoder();
        
        var data = new Uint8Array([0xBA, 0xDD, 0x0B, 0xEE]); // Magic value 0xBADD0BEE (big endian)
        
        for(var i in records) {
            var record = records[i];
            var name = record.name + "." + record.type;
            
            var encoded_name = concatTypedArrays(
                encoder.encode(name),
                new Uint8Array([0])
            );
            
            var encoded_content = concatTypedArrays(
                encoded_name,
                new Uint8Array(await record.data.arrayBuffer())
            );
            
            var length_buffer = new Uint8Array([0xFF, 0xFF]);
            
            encoded_content = concatTypedArrays(length_buffer, encoded_content);
            
            var dv = new DataView(encoded_content.buffer);
            dv.setUint16(0, encoded_content.length, true);
            
            if (data.length + encoded_content.length + 2 > maxSize) {
                console.error("Too much data!");
                throw new Error("Too much data!");
            }
            
            data = concatTypedArrays(data, encoded_content);
        }
        
        data = concatTypedArrays(data, new Uint8Array([0, 0]));
        
        return new Blob([data]);
    }
    
    async __encodeRecord(record) {
        var encoders = this.__getRecordEncoders();
        
        if (record.type in encoders) {
            record = encoders[record.type](record);
        }
        
        return record;
    }
    
    /**
     * Encode the storage from data stored in this class.
     * The second 0xBAD00BEE isn't included.
     *
     * @param   size        max size the storage can take
     *
     * @return  a blob, representing the encoded storage.
     *
     * @throw   Errors      when too much data is passed.
     */
    async encodeStorage(size) {
        var records = Object.assign({}, this.records);
        
        for(var i in this.records) {
            records[i] = await this.__encodeRecord(records[i]);
            
        }
        
        return await this.__assembleStorage(records, size);
    }
    
    async __sliceStorage(blob) {
        var dv = new DataView(await blob.arrayBuffer());
        
        if (dv.getUint32(0x00, false) === 0xBADD0BEE) {
            var offset = 4;
            var records = [];
            
            do {
                var size = dv.getUint16(offset, true);
                
                if (size === 0) break;
                
                var name = this.__readString(dv, offset + 2, size - 2);
                
                var data = blob.slice(offset + 2 + name.size, offset + size);
              //console.log(name,data);
                var record = {
                    name: name.content.split(/\.(?=[^\.]+$)/)[0], // eslint-disable-line no-useless-escape
                    type: name.content.split(/\.(?=[^\.]+$)/)[1], // eslint-disable-line no-useless-escape
                    data: data,
		    size: size,
                }
                
                records.push(record);
                
                offset += size;
                
            } while (size !== 0 && offset < blob.size);
            
            return records;
        } else {
            return {};
        }
    }
    
    __readString(dv, index, maxLen) {
        var out = "";
        var i = 0;
        for(i = 0; i < maxLen || maxLen === 0; i++) {
            var chr = dv.getUint8(index + i);
            
            if (chr === 0) {
                break;
            }
            
            out += String.fromCharCode(chr);
        }
        
        return {
            size: i + 1,
            content: out
        };
    }
    
    async __parsePyRecord(record) {
        var dv = new DataView(await record.data.arrayBuffer());
        
        record.autoImport = dv.getUint8(0) !== 0;
        record.code = this.__readString(dv, 1, record.data.size - 1).content;
        
        //delete record.data;
        
        return record;
    }
    
    __getRecordParsers() {
        return {
            py: this.__parsePyRecord.bind(this)
        };
    }
    
    async __parseRecord(record) {
        var parsers = this.__getRecordParsers();
        
        if (record.type in parsers) {
            record = parsers[record.type](record);
        }
        
        return record;
    }
    
    /**
     * Decode the storage.
     *
     * @param   blob        the encoded storage.
     */
    async parseStorage(blob) {
        var dv = new DataView(await blob.arrayBuffer());
        
        this.magik = dv.getUint32(0x00, false) === 0xBADD0BEE;
    
        this.records = {};
        
        if (this.magik) {
            this.records = await this.__sliceStorage(blob);
            
            for(var i in this.records) {
                this.records[i] = await this.__parseRecord(this.records[i]);
                
                // Throwing away non-python stuff, for convinience.
                // if (this.records[i].type !== 'py') this.records.splice(i, 1);
            }
        }
    }
}

function concatTypedArrays(a, b) {
    // Checks for truthy values on both arrays
    if(!a && !b) throw new Error('Please specify valid arguments for parameters a and b.');  

    // Checks for truthy values or empty arrays on each argument
    // to avoid the unnecessary construction of a new array and
    // the type comparison
    if(!b || b.length === 0) return a;
    if(!a || a.length === 0) return b;

    // Make sure that both typed arrays are of the same type
    if(Object.prototype.toString.call(a) !== Object.prototype.toString.call(b))
        throw new Error('The types of the two arguments passed for parameters a and b do not match.');

    var c = new a.constructor(a.length + b.length);
    c.set(a);
    c.set(b, a.length);

    return c;
}

module.exports = Storage;



/***/ }),
/* 3 */
/***/ (function(module, exports, __webpack_require__) {


const Numworks = __webpack_require__(4);

module.exports = Numworks;



/***/ }),
/* 4 */
/***/ (function(module, exports, __webpack_require__) {


var WebDFU = __webpack_require__(0);
var DFU = WebDFU.DFU;
var DFUse = WebDFU.DFUse;

var Storage = __webpack_require__(2);
var Recovery = __webpack_require__(6);

const AUTOCONNECT_DELAY = 1000;

/**
 * Class handling communication with a Numworks
 * calculator using WebUSB and the WebDFU lib.
 *
 * @author Maxime "M4x1m3" FRIESS
 * @license MIT
 */
class Numworks {
    constructor() {
        this.device = null;
        this.transferSize = 2048;
        this.manifestationTolerant = false;
        this.autoconnectId = null;
    }
    
    /**
     * Get the model of the calculator.
     *
     * @param   exclude_modded  Only include calculator which can be officially purchased from Numworks.
     *                          This includes "0100" and "0110". If a modded Numworks is found, it'll show
     *                          the unmoded version (eg. "0100-8M" becomes "0100").
     *
     * @return  "0110" for an unmodified n0110 (64K internal 8M external).                      "0110" is returned with {exclude_modded}.
     *          "0110-0M" for a modified n0110 (64K internal, no external).                     "????" is returned with {exclude_modded}.
     *          "0110-16M" for a modified n0110 (64K internal, 16M external).                   "0110" is returned with {exclude_modded}.
     *          "0100" for unmodified n0100 (1M internal, no external).                         "0100" is returned with {exclude_modded}.
     *          "0100-8M"  for a "Numworks++" with 8M external (1M internal, 8M external).      "0100" is returned with {exclude_modded}.
     *          "0100-16M" for a "Numworks++" with 16M external (1M internal, 16M external).    "0100" is returned with {exclude_modded}.
     *
     *          Other flash sizes don't exist for the packaging the Numworks (SOIC-8) uses, so it's safe to assume
     *          we'll only encounter 0M, 8M and 16M versions.
     *
     *          "????" if can't be determined (maybe the user plugged a DFU capable device which isn't a Numworks).
     */
    getModel(exclude_modded = true) {
        var internal_size = 0;
        var external_size = 0;
        
        for (let i = 0; i < this.device.memoryInfo.segments.length; i++) {
            
            if (this.device.memoryInfo.segments[i].start >= 0x08000000 && this.device.memoryInfo.segments[i].start <= 0x080FFFFF) {
                internal_size += this.device.memoryInfo.segments[i].end - this.device.memoryInfo.segments[i].start;
            }
            
            if (this.device.memoryInfo.segments[i].start >= 0x90000000 && this.device.memoryInfo.segments[i].start <= 0x9FFFFFFF) {
                external_size += this.device.memoryInfo.segments[i].end - this.device.memoryInfo.segments[i].start;
            }
        }
        
        if (internal_size === 0x10000) {
            if (external_size === 0) {
                return (exclude_modded ? "????" : "0110-0M");
            } else if (external_size === 0x800000) {
                return "0110";
            } else if (external_size === 0x1000000) {
                return (exclude_modded ? "0110" : "0110-16M");
            } else {
                return "????";
            }
        } else if (internal_size === 0x100000) {
            if (external_size === 0) {
                return "0100";
            } else if (external_size === 0x800000) {
                return (exclude_modded ? "0100" : "0100-8M");
            } else if (external_size === 0x1000000) {
                return (exclude_modded ? "0100" : "0100-16M");
            } else {
                return "????";
            }
        } else {
            return "????";
        }
    }
    
    /**
     * Flash buffer to internal flash.
     *
     * @param   buffer      ArrayBuffer to flash.
     */
    async flashInternal(buffer) {
      this.device.startAddress = 0x08000000;
      this.device.logProgress(0,'Erase and flash internal');
      await this.device.do_download(this.transferSize, buffer, true);
    }
    
    /**
     * Flash buffer to external flash.
     *
     * @param   buffer      ArrayBuffer to flash.
     */
    async flashExternal(buffer) {
      this.device.startAddress = 0x90000000;
      this.device.logProgress(0,'Erase and flash external');
      await this.device.do_download(this.transferSize, buffer, false);
    }
    
    async __getDFUDescriptorProperties(device) {
        // Attempt to read the DFU functional descriptor
        // TODO: read the selected configuration's descriptor
        return device.readConfigurationDescriptor(0).then(
            data => {
                let configDesc = DFU.parseConfigurationDescriptor(data);
                let funcDesc = null;
                let configValue = device.settings.configuration.configurationValue;
                if (configDesc.bConfigurationValue === configValue) {
                    for (let desc of configDesc.descriptors) {
                        if (desc.bDescriptorType === 0x21 && desc.hasOwnProperty("bcdDFUVersion")) {
                            funcDesc = desc;
                            break;
                        }
                    }
                }

                if (funcDesc) {
                    return {
                        WillDetach:            ((funcDesc.bmAttributes & 0x08) !== 0),
                        ManifestationTolerant: ((funcDesc.bmAttributes & 0x04) !== 0),
                        CanUpload:             ((funcDesc.bmAttributes & 0x02) !== 0),
                        CanDnload:             ((funcDesc.bmAttributes & 0x01) !== 0),
                        TransferSize:          funcDesc.wTransferSize,
                        DetachTimeOut:         funcDesc.wDetachTimeOut,
                        DFUVersion:            funcDesc.bcdDFUVersion
                    };
                } else {
                    return {};
                }
            },
            error => {}
        );
    }
    
    /**
     * Detect a numworks calculator.
     *
     * @param   successCallback     Callback in case of success.
     * @param   errorCallback       Callback in case of error.
     */
    async detect(successCallback, errorCallback) {
        var _this = this;
        navigator.usb.requestDevice({ 'filters': [{'vendorId': 0x0483, 'productId': 0xa291}]}).then(
            async selectedDevice => {
                let interfaces = DFU.findDeviceDfuInterfaces(selectedDevice);
                await _this.__fixInterfaceNames(selectedDevice, interfaces);
                _this.device = await _this.__connect(new DFU.Device(selectedDevice, interfaces[0]));
                
                successCallback();
            }
        ).catch(error => {
            errorCallback(error);
        });
    }
    
    /**
     * Connect to a WebDFU device.
     *
     * @param   device      The WebUSB device to connect to.
     */
    async __connect(device) {
        try {
            await device.open();
        } catch (error) {
            // this.installInstance.calculatorError(true, error);
            throw error;
        }

        // Attempt to parse the DFU functional descriptor
        let desc = {};
        try {
            desc = await this.__getDFUDescriptorProperties(device);
        } catch (error) {
            // this.installInstance.calculatorError(true, error);
            throw error;
        }

        if (desc && Object.keys(desc).length > 0) {
            device.properties = desc;
            this.transferSize = desc.TransferSize;
            if (desc.CanDnload) {
                this.manifestationTolerant = desc.ManifestationTolerant;
            }

            if ((desc.DFUVersion === 0x100 || desc.DFUVersion === 0x011a) && device.settings.alternate.interfaceProtocol === 0x02) {
                device = new DFUse.Device(device.device_, device.settings);
                if (device.memoryInfo) {
                    // We have to add RAM manually, because the device doesn't expose that normally
                    device.memoryInfo.segments.unshift({
                        start: 0x20000000,
                        sectorSize: 1024,
                        end: 0x20040000,
                        readable: true,
                        erasable: false,
                        writable: true
                    });
                }
            }
        }

        // Bind logging methods
        device.logDebug = console.log;
        device.logInfo = console.info;
        device.logWarning = console.warn;
        device.logError = console.error;
        //device.logProgress = console.log;
        
        return device;
    }
    
    __readFString(dv, index, len) {
        var out = "";
        for(var i = 0; i < len; i++) {
            var chr = dv.getUint8(index + i);
            
            if (chr === 0) {
                break;
            }
            
            out += String.fromCharCode(chr);
        }
        
        return out;
    }
    
  __parsePlatformInfo(array,internalflash) {
    var dv = new DataView(array);
    var data = {};
    
    //data["magik"] = dv.getUint32(0x00, false) === 0xF00DC0DE;
    
    data["magik"] = dv.getUint32(0x00, false) === 0xF00DC0DE;
    if (internalflash){
      if (data["magik"]) {
        data["oldplatform"] = !(dv.getUint32(0x1C, false) === 0xF00DC0DE);
        
        data["omega"] = {};
        
        if (data["oldplatform"]) {
          data["omega"]["installed"] = dv.getUint32(0x1C + 8, false) === 0xF00DC0DE || dv.getUint32(0x1C + 16, false) === 0xDEADBEEF || dv.getUint32(0x1C + 32, false) === 0xDEADBEEF;
          if (data["omega"]["installed"]) {
            data["omega"]["version"] = this.__readFString(dv, 0x0C, 16);
            
            data["omega"]["user"] = "";
            
          }
          
          data["version"] = this.__readFString(dv, 0x04, 8);
          var offset = 0;
          if (dv.getUint32(0x1C + 8, false) === 0xF00DC0DE) {
            offset = 8;
          } else if (dv.getUint32(0x1C + 16, false) === 0xF00DC0DE) {
            offset = 16;
          } else if (dv.getUint32(0x1C + 32, false) === 0xF00DC0DE) {
            offset = 32;
          }
          
          data["commit"] = this.__readFString(dv, 0x0C + offset, 8);
          data["storage"] = {};
          data["storage"]["address"] = dv.getUint32(0x14 + offset, true);
          data["storage"]["size"] = dv.getUint32(0x18 + offset, true);
        } else {
          data["omega"]["installed"] = dv.getUint32(0x20, false) === 0xDEADBEEF && dv.getUint32(0x44, false) === 0xDEADBEEF;
          if (data["omega"]["installed"]) {
            data["omega"]["version"] = this.__readFString(dv, 0x24, 16);
            data["omega"]["user"] = this.__readFString(dv, 0x34, 16);
          }
	  
          data["version"] = this.__readFString(dv, 0x04, 8);
          data["commit"] = this.__readFString(dv, 0x0C, 8);
          data["storage"] = {};
          data["storage"]["address"] = dv.getUint32(0x14, true);
          data["storage"]["size"] = dv.getUint32(0x18, true);
        }
      } else {
        data["omega"] = false;
      }
      
      return data;
    } // end internalflash
    data["storage"] = {};
    data["storage"]["address"] = dv.getUint32(0xc, true);
    data["storage"]["size"] = dv.getUint32(0x10, true);
    // console.log('external flash platinfo',data["storage"]["address"],data["storage"]["size"]);
    return data;
  }
  
    
    /**
     * Get the platforminfo section of the calculator.
     *
     * @return  an object representing the platforminfo.
     */
  async getPlatformInfo() {
    let tab=[0x90000000,0x90180000,0x90400000],res;
    for (let j=0;j<tab.length;j++){
      this.device.startAddress = tab[j];
      const blobk = await this.device.do_upload(this.transferSize, 0x20);
      let dv=new DataView(await blobk.arrayBuffer()),res;
      console.log('getplatfrominfo 1',tab[j],dv.getUint32(0x0,false),dv.getUint32(0x4,false),dv.getUint32(0x8,false),dv.getUint32(0xc,false));
      if (dv.getUint32(0x8,false)===0xF00DC0DE){ // f0 0d c0 de
	this.device.startAddress += 0x10000;
	const blob = await this.device.do_upload(this.transferSize, 0x48);
	const blobb = await blob.arrayBuffer();
	let dvb=new DataView(blobb);
	console.log('getplatfrominfo 2',tab[j],dv.getUint32(0x0,false),dv.getUint32(0x4,false),dv.getUint32(0x8,false),dv.getUint32(0xc,false));
	this.device.startAddress=dvb.getUint32(0xc, true)
	const blob2 = await this.device.do_upload(this.transferSize, 0x4);
	const blobb2 = await blob2.arrayBuffer();
	let dv2b=new DataView(blobb2);
	console.log(dv2b.getUint32(0x00, false));
	if (dv2b.getUint32(0x00, false) === 0xbadd0bee){
          res=this.__parsePlatformInfo(blobb,false);
	  res["bootloader"]=1;
	  return res;
	}
      }
    }
    this.device.startAddress = 0x080001c4;
    const blob = await this.device.do_upload(this.transferSize, 0x48);
    res=this.__parsePlatformInfo(await blob.arrayBuffer(),true);
    res["bootloader"]=0;
    return res;
  }
    
    async get_internal_flash() {
      this.device.startAddress = 0x08000000;
      this.device.logProgress(0,'Get internal');
      const blob = await this.device.do_upload(this.transferSize, 0x8000);
      return await blob.arrayBuffer();
    }
    
  async get_external_flash(slot) {
    if (slot==2)
      this.device.startAddress = 0x90180000;
    else if (slot==3)
      this.device.startAddress = 0x90400000;
    else
      this.device.startAddress = 0x90000000;
    this.device.logProgress(0,'Get external');
    const blob = await this.device.do_upload(this.transferSize, slot==2?0x80000:0x130000);
    return await blob.arrayBuffer();
  }

    async rw_check(addr,n){
      this.device.startAddress = addr;
      let data=new Uint8Array(n);
      let i;
      for (i=0;i<n;++i){
	let j=Math.floor(Math.random()*n);
	data[j]=Math.floor(Math.random()*0xff);
      }
      this.device.logProgress(0,'Read/Write check');
      await this.device.do_download(this.transferSize, data, false);
      const blob = await this.device.do_upload(this.transferSize, n);
      let data1= new Uint8Array(await blob.arrayBuffer());
      for (i=0;i<n;++i){
	if (data[i]!=data1[i])
	  break;
      }
      return i==n;
    }
  
    async get_apps() {
      this.device.startAddress = 0x90200000;
      this.device.logProgress(0,'Get apps');
      const blob = await this.device.do_upload(this.transferSize, 0x600000);
      return await blob.arrayBuffer();
    }
    
    async __autoConnectDevice(device) {
        let interfaces = DFU.findDeviceDfuInterfaces(device.device_);
        await this.__fixInterfaceNames(device.device_, interfaces);
        device = await this.__connect(new DFU.Device(device.device_, interfaces[0]));
        return device;
    }
    
    /**
     * Autoconnect a numworks calculator
     *
     * @param   serial      Serial number. If ommited, any will work.
     */
    autoConnect(callback, serial) {
        var _this = this;
        var vid = 0x0483, pid = 0xa291;
      console.log('autoconnect',DFU);
        DFU.findAllDfuInterfaces().then(async dfu_devices => {
          console.log('autoconnect',dfu_devices);
            let matching_devices = _this.__findMatchingDevices(vid, pid, serial, dfu_devices);
            
            if (matching_devices.length !== 0) {
                this.stopAutoConnect();
                
                this.device = await this.__autoConnectDevice(matching_devices[0]);
                
                await callback();
            }
        });
        
        this.autoconnectId = setTimeout(this.autoConnect.bind(this, callback, serial), AUTOCONNECT_DELAY);
    }
    
    /**
     * Stop autoconnection.
     */
    stopAutoConnect() {
        if (this.autoconnectId === null) return;
        
        clearTimeout(this.autoconnectId);
        
        this.autoconnectId = null;
    }
    
    async __fixInterfaceNames(device_, interfaces) {
        // Check if any interface names were not read correctly
        if (interfaces.some(intf => (intf.name === null))) {
            // Manually retrieve the interface name string descriptors
            let tempDevice = new DFU.Device(device_, interfaces[0]);
            await tempDevice.device_.open();
            let mapping = await tempDevice.readInterfaceNames();
            await tempDevice.close();

            for (let intf of interfaces) {
                if (intf.name === null) {
                    let configIndex = intf.configuration.configurationValue;
                    let intfNumber = intf["interface"].interfaceNumber;
                    let alt = intf.alternate.alternateSetting;
                    intf.name = mapping[configIndex][intfNumber][alt];
                }
            }
        }
    }
    
    __findMatchingDevices(vid, pid, serial, dfu_devices) {
        let matching_devices = [];
        for (let dfu_device of dfu_devices) {
            if (serial) {
                if (dfu_device.device_.serialNumber === serial) {
                    matching_devices.push(dfu_device);
                }
            } else {
                if (
                    (!pid && vid > 0 && dfu_device.device_.vendorId  === vid) ||
                    (!vid && pid > 0 && dfu_device.device_.productId === pid) ||
                    (vid > 0 && pid > 0 && dfu_device.device_.vendorId  === vid && dfu_device.device_.productId === pid)
                   )
                {
                    matching_devices.push(dfu_device);
                }
            }
        }
        
        return matching_devices;
    }
    
    /**
     * Get storage from the calculator.
     *
     * @param   address     Storage address
     * @param   size        Storage size.
     *
     * @return  The sotrage, as a Blob.
     */
    async __retreiveStorage(address, size) {
        this.device.startAddress = address;
        return await this.device.do_upload(this.transferSize, size + 8);
    }
    
    /**
     * Flash storage to the calculator.
     *
     * @param   address     Storage address
     * @param   data        Storage data.
     */
    async __flashStorage(address, data) {
        this.device.startAddress = address;
        await this.device.do_download(this.transferSize, data, false);
    }
    
    /**
     * Install new storage in calculator
     *
     * @param   storage     Storage class, representing the storage.
     * @param   callback    Callback to be called when done.
     *
     * @throw   Error       If storage is too big.
     */
    async installStorage(storage, callback) {
        let pinfo = await this.getPlatformInfo();
        
        let storage_blob = await storage.encodeStorage(pinfo["storage"]["size"]);
        await this.__flashStorage(pinfo["storage"]["address"], await storage_blob.arrayBuffer());
        
        callback();
    }
    
    /**
     * Get and parse storage on the calculator.
     *
     * @return  Storage class describing the storage of the calculator.
     */
    async backupStorage() {
        let pinfo = await this.getPlatformInfo();
        
        let storage_blob = await this.__retreiveStorage(pinfo["storage"]["address"], pinfo["storage"]["size"]);
        
        let storage = new Numworks.Storage();
        
        await storage.parseStorage(storage_blob);
        
        return storage;
    }
    
    onUnexpectedDisconnect(event, callback) {
        if (this.device !== null && this.device.device_ !== null) {
            if (this.device.device_ === event.device) {
                this.device.disconnected = true;
                callback(event);
                this.device = null;
            }
        }
    }
}

Numworks.Recovery = Recovery;
Numworks.Storage = Storage;

module.exports = Numworks;



/***/ }),
/* 5 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


const DFU = __webpack_require__(1);

class DFUse extends DFU {
    static get GET_COMMANDS() { return 0x00; }
    static get SET_ADDRESS() { return 0x21; }
    static get ERASE_SECTOR() { return 0x41; }
    
    static parseMemoryDescriptor(desc) {
        const nameEndIndex = desc.indexOf("/");
        if (!desc.startsWith("@") || nameEndIndex === -1) {
            throw new Error(`Not a DfuSe memory descriptor: "${desc}"`);
        }

        const name = desc.substring(1, nameEndIndex).trim();
        const segmentString = desc.substring(nameEndIndex);

        let segments = [];

        const sectorMultipliers = {
            ' ': 1,
            'B': 1,
            'K': 1024,
            'M': 1048576
        };

        let contiguousSegmentRegex = /\/\s*(0x[0-9a-fA-F]{1,8})\s*\/(\s*[0-9]+\s*\*\s*[0-9]+\s?[ BKM]\s*[abcdefg]\s*,?\s*)+/g;
        let contiguousSegmentMatch;
        while ((contiguousSegmentMatch = contiguousSegmentRegex.exec(segmentString))) {
            let segmentRegex = /([0-9]+)\s*\*\s*([0-9]+)\s?([ BKM])\s*([abcdefg])\s*,?\s*/g;
            let startAddress = parseInt(contiguousSegmentMatch[1], 16);
            let segmentMatch;
            while ((segmentMatch = segmentRegex.exec(contiguousSegmentMatch[0]))) {
                let segment = {}
                let sectorCount = parseInt(segmentMatch[1], 10);
                let sectorSize = parseInt(segmentMatch[2]) * sectorMultipliers[segmentMatch[3]];
                let properties = segmentMatch[4].charCodeAt(0) - 'a'.charCodeAt(0) + 1;
                segment.start = startAddress;
                segment.sectorSize = sectorSize;
                segment.end = startAddress + sectorSize * sectorCount;
                segment.readable = (properties & 0x1) !== 0;
                segment.erasable = (properties & 0x2) !== 0;
                segment.writable = (properties & 0x4) !== 0;
                segments.push(segment);

                startAddress += sectorSize * sectorCount;
            }
        }

        return {"name": name, "segments": segments};
    }
}

DFUse.Device = class extends DFU.Device {
    constructor(device, settings) {
        super(device, settings);
        this.memoryInfo = null;
        this.startAddress = NaN;
        if (settings.name) {
            this.memoryInfo = DFUse.parseMemoryDescriptor(settings.name);
        }
    }
    
    async dfuseCommand(command, param, len) {
        if (typeof param === 'undefined' && typeof len === 'undefined') {
            param = 0x00;
            len = 1;
        }

        const commandNames = {
            0x00: "GET_COMMANDS",
            0x21: "SET_ADDRESS",
            0x41: "ERASE_SECTOR"
        };

        let payload = new ArrayBuffer(len + 1);
        let view = new DataView(payload);
        view.setUint8(0, command);
        if (len === 1) {
            view.setUint8(1, param);
        } else if (len === 4) {
            view.setUint32(1, param, true);
        } else {
            throw new Error("Don't know how to handle data of len " + len);
        }

        try {
            await this.download(payload, 0);
        } catch (error) {
            throw new Error("Error during special DfuSe command " + commandNames[command] + ":" + error);
        }

        let status = await this.poll_until(state => (state !== DFU.dfuDNBUSY));
        if (status.status !== DFU.STATUS_OK) {
            throw new Error("Special DfuSe command " + commandNames[command] + " failed");
        }
    }
    
    getSegment(addr) {
        if (!this.memoryInfo || ! this.memoryInfo.segments) {
            throw new Error("No memory map information available");
        }
        for (let segment of this.memoryInfo.segments) {
            if (segment.start <= addr && addr < segment.end) {
	      // console.log('getSegment',segment,addr);
                return segment;
            }
        }

        return null;
    }
    
    getSectorStart(addr, segment) {
        if (typeof segment === 'undefined') {
            segment = this.getSegment(addr);
        }

        if (!segment) {
            throw new Error(`Address ${addr.toString(16)} outside of memory map`);
        }

        const sectorIndex = Math.floor((addr - segment.start)/segment.sectorSize);
        return segment.start + sectorIndex * segment.sectorSize;
    }

  getSectorEnd(addr, segment) {
    //console.log(addr,segment);console.trace();
        if (typeof segment === 'undefined') {
            segment = this.getSegment(addr);
        }

        if (!segment) {
            throw new Error(`Address ${addr.toString(16)} outside of memory map`);
        }

        const sectorIndex = Math.floor((addr - segment.start)/segment.sectorSize);
        return segment.start + (sectorIndex + 1) * segment.sectorSize;
    }

    getFirstWritableSegment() {
        if (!this.memoryInfo || ! this.memoryInfo.segments) {
            throw new Error("No memory map information available");
        }

        for (let segment of this.memoryInfo.segments) {
            if (segment.writable) {
                return segment;
            }
        }

        return null;
    }

    getMaxReadSize(startAddr) {
        if (!this.memoryInfo || ! this.memoryInfo.segments) {
            throw new Error("No memory map information available");
        }

        let numBytes = 0;
        for (let segment of this.memoryInfo.segments) {
            if (segment.start <= startAddr && startAddr < segment.end) {
                // Found the first segment the read starts in
                if (segment.readable) {
                    numBytes += segment.end - startAddr;
                } else {
                    return 0;
                }
            } else if (segment.start === startAddr + numBytes) {
                // Include a contiguous segment
                if (segment.readable) {
                    numBytes += (segment.end - segment.start);
                } else {
                    break;
                }
            }
        }

        return numBytes;
    };

  async erase(startAddr, length) {
    let segment = this.getSegment(startAddr);
    console.log('erase',startAddr,length,segment);
    if (!segment.erasable){
      if (segment.writable) return 1;
      return 0;
    }
    let addr = this.getSectorStart(startAddr, segment);
    const endAddr = this.getSectorEnd(startAddr + length - 1);
    
    let bytesErased = 0;
    const bytesToErase = endAddr - addr;
    if (bytesToErase > 0) {
      this.logProgress(bytesErased, bytesToErase);
    }
    
    while (addr < endAddr) {
      if (segment.end <= addr) {
        segment = this.getSegment(addr);
      }
      if (!segment.erasable) {
        // Skip over the non-erasable section
        bytesErased = Math.min(bytesErased + segment.end - addr, bytesToErase);
        addr = segment.end;
        this.logProgress(bytesErased, bytesToErase);
        continue;
      }
      const sectorIndex = Math.floor((addr - segment.start)/segment.sectorSize);
      const sectorAddr = segment.start + sectorIndex * segment.sectorSize;
      this.logDebug(`Erasing ${segment.sectorSize}B at 0x${sectorAddr.toString(16)}`);
      await this.dfuseCommand(DFUse.ERASE_SECTOR, sectorAddr, 4);
      addr = sectorAddr + segment.sectorSize;
      bytesErased += segment.sectorSize;
      this.logProgress(bytesErased, bytesToErase);
    }
    return 1;
  };

    async do_download(xfer_size, data, manifestationTolerant) {
      if (!this.memoryInfo || ! this.memoryInfo.segments) {
        throw new Error("No memory map available");
      }
      
      this.logInfo("Erasing DFU device memory");
      
      let bytes_sent = 0;
      let expected_size = data.byteLength;
      
      let startAddress = this.startAddress;
      if (isNaN(startAddress)) {
        startAddress = this.memoryInfo.segments[0].start;
        this.logWarning("Using inferred start address 0x" + startAddress.toString(16));
      } else if (this.getSegment(startAddress) === null) {
        this.logError(`Start address 0x${startAddress.toString(16)} outside of memory map bounds`);
	return 0;
      }
      let savelegend='';
      if (typeof UI!=='undefined' && typeof UI.nws_progresslegend!=='undefined'){
	savelegend=UI.nws_progresslegend.innerHTML;
	this.logProgress(0,'Erase '+savelegend);
	//console.log('do_download before erase',savelegend);
      }
      let res=0;
      try {
	res=await this.erase(startAddress, expected_size);
      } catch(err) { res=0;}
      if (res==0) return 0;
      if (savelegend!=''){
	//console.log('do_download after erase',savelegend);
	this.logProgress(0,'Flash '+savelegend);
      }
      this.logInfo("Copying data from browser to DFU device");
      
      let address = startAddress;
      while (bytes_sent < expected_size) {
        const bytes_left = expected_size - bytes_sent;
        const chunk_size = Math.min(bytes_left, xfer_size);
	
        let bytes_written = 0;
        let dfu_status;
        try {
          await this.dfuseCommand(DFUse.SET_ADDRESS, address, 4);
          this.logDebug(`Set address to 0x${address.toString(16)}`);
          bytes_written = await this.download(data.slice(bytes_sent, bytes_sent+chunk_size), 2);
          this.logDebug("Sent " + bytes_written + " bytes");
          dfu_status = await this.poll_until_idle(DFU.dfuDNLOAD_IDLE);
          address += chunk_size;
        } catch (error) {
          throw new Error("Error during DfuSe download: " + error);
        }
	
        if (dfu_status.status !== DFU.STATUS_OK) {
          throw new Error(`DFU DOWNLOAD failed state=${dfu_status.state}, status=${dfu_status.status}`);
        }
	
        this.logDebug("Wrote " + bytes_written + " bytes");
        bytes_sent += bytes_written;
	
        this.logProgress(bytes_sent, expected_size);
      }
      this.logInfo(`Wrote ${bytes_sent} bytes`);
      
      if(manifestationTolerant) {
        this.logInfo("Manifesting new firmware");
        try {
          await this.dfuseCommand(DFUse.SET_ADDRESS, startAddress, 4);
          await this.download(new ArrayBuffer(), 2);
        } catch (error) {
          throw new Error("Error during DfuSe manifestation: " + error);
        }
	
        try {
          await this.poll_until(state => (state === DFU.dfuMANIFEST));
        } catch (error) {
          this.logError(error);
        }
      }
      return 1;
    }
  
    async do_upload(xfer_size, max_size) {
        let startAddress = this.startAddress;
        if (isNaN(startAddress)) {
            startAddress = this.memoryInfo.segments[0].start;
            this.logWarning("Using inferred start address 0x" + startAddress.toString(16));
        } else if (this.getSegment(startAddress) === null) {
            this.logWarning(`Start address 0x${startAddress.toString(16)} outside of memory map bounds`);
        }

        this.logInfo(`Reading up to 0x${max_size.toString(16)} bytes starting at 0x${startAddress.toString(16)}`);
        let state = await this.getState();
        if (state !== DFU.dfuIDLE) {
            await this.abortToIdle();
        }
        await this.dfuseCommand(DFUse.SET_ADDRESS, startAddress, 4);
        await this.abortToIdle();

        // DfuSe encodes the read address based on the transfer size,
        // the block number - 2, and the SET_ADDRESS pointer.
        return await super.do_upload(xfer_size, max_size, 2);
    }
};

module.exports = DFUse;




/***/ }),
/* 6 */
/***/ (function(module, exports, __webpack_require__) {


var WebDFU = __webpack_require__(0);
var DFU = WebDFU.DFU;
var DFUse = WebDFU.DFUse;

var Storage = __webpack_require__(2);

const AUTOCONNECT_DELAY = 1000;

/**
 * Class handling communication with a Numworks
 * calculator in Recovery Mode using WebUSB and the WebDFU lib.
 *
 * @author Maxime "M4x1m3" FRIESS
 * @license MIT
 */
class Recovery {
    constructor() {
        this.device = null;
        this.transferSize = 2048;
        this.manifestationTolerant = false;
        this.autoconnectId = null;
    }
    
    /**
     * Get approximated the model of the calculator.
      
     * This just checks the size of the internal size, because that's everything the STM32 bootloader
     * exposes.
     *
     * @note    The check for the N0110 **WILL** break if a new model happens to actually have 512K internal.
     *          We have to ckeck for 512K because every STM32F73x bootloaders advertize 512K regardless of
     *          the actual capacity of the internal flah.
     *          TODO: Find a better way to detect the model (Numworks' private API ?)
     *
     * @return  "0110" for an unmodified n0110 (512K advertized internal).
     *          "0100" for unmodified n0100 (1M internal).
     *          "????" for something unknown (Other internal sizes).
     */
    getModel(exclude_modded = true) {
        var internal_size = 0;
        
        for (let i = 0; i < this.device.memoryInfo.segments.length; i++) {
            if (this.device.memoryInfo.segments[i].start >= 0x08000000 && this.device.memoryInfo.segments[i].start <= 0x080FFFFF) {
                internal_size += this.device.memoryInfo.segments[i].end - this.device.memoryInfo.segments[i].start;
            }
        }
        
        if (internal_size === 0x80000) {
            return "0110";
        } else if (internal_size === 0x100000) {
            return "0100";
        } else {
            return "????";
        }
    }
    
    /**
     * Flash buffer to recovery location, in RAM.
     *
     * @param   buffer      ArrayBuffer to flash.
     */
    async flashRecovery(buffer) {
        this.device.startAddress = 0x20030000;
        // This is needed because the STM32F73x bootloader starts in dfuERROR status,
        // for a weird reason that I spend hours figuring out, but didn't find.
        // Better to not think about it.
        await this.device.clearStatus();
        await this.device.do_download(this.transferSize, buffer, true);
    }
    
    async __getDFUDescriptorProperties(device) {
        // Attempt to read the DFU functional descriptor
        // TODO: read the selected configuration's descriptor
        return device.readConfigurationDescriptor(0).then(
            data => {
                let configDesc = DFU.parseConfigurationDescriptor(data);
                let funcDesc = null;
                let configValue = device.settings.configuration.configurationValue;
                if (configDesc.bConfigurationValue === configValue) {
                    for (let desc of configDesc.descriptors) {
                        if (desc.bDescriptorType === 0x21 && desc.hasOwnProperty("bcdDFUVersion")) {
                            funcDesc = desc;
                            break;
                        }
                    }
                }

                if (funcDesc) {
                    return {
                        WillDetach:            ((funcDesc.bmAttributes & 0x08) !== 0),
                        ManifestationTolerant: ((funcDesc.bmAttributes & 0x04) !== 0),
                        CanUpload:             ((funcDesc.bmAttributes & 0x02) !== 0),
                        CanDnload:             ((funcDesc.bmAttributes & 0x01) !== 0),
                        TransferSize:          funcDesc.wTransferSize,
                        DetachTimeOut:         funcDesc.wDetachTimeOut,
                        DFUVersion:            funcDesc.bcdDFUVersion
                    };
                } else {
                    return {};
                }
            },
            error => {}
        );
    }
    
    /**
     * Detect a numworks calculator.
     *
     * @param   successCallback     Callback in case of success.
     * @param   errorCallback       Callback in case of error.
     */
    async detect(successCallback, errorCallback) {
        var _this = this;
        navigator.usb.requestDevice({ 'filters': [{'vendorId': 0x0483, 'productId': 0xdf11}]}).then(
            async selectedDevice => {
                let interfaces = DFU.findDeviceDfuInterfaces(selectedDevice);
                await _this.__fixInterfaceNames(selectedDevice, interfaces);
                _this.device = await _this.__connect(new DFU.Device(selectedDevice, interfaces[0]));
                
                successCallback();
            }
        ).catch(error => {
            errorCallback(error);
        });
    }
    
    /**
     * Connect to a WebDFU device.
     *
     * @param   device      The WebUSB device to connect to.
     */
    async __connect(device) {
        try {
            await device.open();
        } catch (error) {
            // this.installInstance.calculatorError(true, error);
            throw error;
        }

        // Attempt to parse the DFU functional descriptor
        let desc = {};
        try {
            desc = await this.__getDFUDescriptorProperties(device);
        } catch (error) {
            // this.installInstance.calculatorError(true, error);
            throw error;
        }

        if (desc && Object.keys(desc).length > 0) {
            device.properties = desc;
            this.transferSize = desc.TransferSize;
            if (desc.CanDnload) {
                this.manifestationTolerant = desc.ManifestationTolerant;
            }

            if ((desc.DFUVersion === 0x100 || desc.DFUVersion === 0x011a) && device.settings.alternate.interfaceProtocol === 0x02) {
                device = new DFUse.Device(device.device_, device.settings);
                if (device.memoryInfo) {
                    // We have to add RAM manually, because the device doesn't expose that normally
                    device.memoryInfo.segments.unshift({
                        start: 0x20000000,
                        sectorSize: 1024,
                        end: 0x20040000,
                        readable: true,
                        erasable: false,
                        writable: true
                    });
                }
            }
        }

        // Bind logging methods
        device.logDebug = console.log;
        device.logInfo = console.info;
        device.logWarning = console.warn;
        device.logError = console.error;
        device.logProgress = console.log;
        
        return device;
    }
    
    async __autoConnectDevice(device) {
        let interfaces = DFU.findDeviceDfuInterfaces(device.device_);
        await this.__fixInterfaceNames(device.device_, interfaces);
        device = await this.__connect(new DFU.Device(device.device_, interfaces[0]));
        return device;
    }
    
    /**
     * Autoconnect a numworks calculator
     *
     * @param   serial      Serial number. If ommited, any will work.
     */
    async autoConnect(callback, serial) {
        var _this = this;
        var vid = 0x0483, pid = 0xdf11;
      console.log('rescue autoconnect',DFU);
      var chk=await navigator.usb.getDevices();
      if (chk.length==0)
      {
	await navigator.usb.requestDevice({ filters: [{ vendorId: 0x0483 }] });	
	//return -1;
      }
      DFU.findAllDfuInterfaces().then(async dfu_devices => {
          console.log('rescue autoconnect',dfu_devices);
            let matching_devices = _this.__findMatchingDevices(vid, pid, serial, dfu_devices);
            if (matching_devices.length !== 0) {
                this.stopAutoConnect();
                
                this.device = await this.__autoConnectDevice(matching_devices[0]);
                
                await callback();
            }
        });
        
      this.autoconnectId = setTimeout(this.autoConnect.bind(this, callback, serial), AUTOCONNECT_DELAY);
      return 0;
    }
    
    /**
     * Stop autoconnection.
     */
  stopAutoConnect() {
    console.log('rescue stop auto',this.autoconnectId);
        if (this.autoconnectId === null) return;
        
        clearTimeout(this.autoconnectId);
        
        this.autoconnectId = null;
    }
    
    async __fixInterfaceNames(device_, interfaces) {
        // Check if any interface names were not read correctly
        if (interfaces.some(intf => (intf.name === null))) {
            // Manually retrieve the interface name string descriptors
            let tempDevice = new DFU.Device(device_, interfaces[0]);
            await tempDevice.device_.open();
            let mapping = await tempDevice.readInterfaceNames();
            await tempDevice.close();

            for (let intf of interfaces) {
                if (intf.name === null) {
                    let configIndex = intf.configuration.configurationValue;
                    let intfNumber = intf["interface"].interfaceNumber;
                    let alt = intf.alternate.alternateSetting;
                    intf.name = mapping[configIndex][intfNumber][alt];
                }
            }
        }
    }
    
    __findMatchingDevices(vid, pid, serial, dfu_devices) {
        let matching_devices = [];
        for (let dfu_device of dfu_devices) {
            if (serial) {
                if (dfu_device.device_.serialNumber === serial) {
                    matching_devices.push(dfu_device);
                }
            } else {
                if (
                    (!pid && vid > 0 && dfu_device.device_.vendorId  === vid) ||
                    (!vid && pid > 0 && dfu_device.device_.productId === pid) ||
                    (vid > 0 && pid > 0 && dfu_device.device_.vendorId  === vid && dfu_device.device_.productId === pid)
                   )
                {
                    matching_devices.push(dfu_device);
                }
            }
        }
        
        return matching_devices;
    }
    
    /**
     * Get storage from the calculator.
     *
     * @param   address     Storage address
     * @param   size        Storage size.
     *
     * @return  The sotrage, as a Blob.
     */
    async __retreiveStorage(address, size) {
        this.device.startAddress = address;
        return await this.device.do_upload(this.transferSize, size + 8);
    }
    
    /**
     * Flash storage to the calculator.
     *
     * @param   address     Storage address
     * @param   data        Storage data.
     */
    async __flashStorage(address, data) {
        this.device.startAddress = address;
        await this.device.do_download(this.transferSize, data, false);
    }
    
    onUnexpectedDisconnect(event, callback) {
        if (this.device !== null && this.device.device_ !== null) {
            if (this.device.device_ === event.device) {
                this.device.disconnected = true;
                callback(event);
                this.device = null;
            }
        }
    }
}

module.exports = Recovery;



/***/ })
/******/ ]);
