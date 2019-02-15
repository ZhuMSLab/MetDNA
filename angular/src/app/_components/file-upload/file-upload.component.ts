// created by Xin Xiong<onlybelter@outlook.com>, https://github.com/OnlyBelter

import { Component, OnInit, EventEmitter, Output, ViewChild, ElementRef } from '@angular/core';
// import { Headers, RequestOptions } from "@angular/http";
import {FormGroup, FormBuilder, Validators} from '@angular/forms';
import {DomSanitizer} from '@angular/platform-browser';

import { UploadOutput, UploadFile, humanizeBytes, UploadInput, UploaderOptions} from 'ngx-uploader';
import { UserService, ProjectService, FileService } from '../../_service/index';
import { sample } from 'rxjs/operator/sample';
import { element } from 'protractor';
import { COMPOSITION_BUFFER_MODE } from '@angular/forms/src/directives/default_value_accessor';
declare var $: any;

const ERROR_CODE = {0: 'valid',
                    1: '[ERROR 1] - ' + 'Missing values are found in the uploaded file',
                    2: '[ERROR 2] - ' + 'Column name of the first column in MS1 peak table is not "name"',
                    3: '[ERROR 3] - ' + 'Column name of the second column in MS1 peak table is not "mz"',
                    4: '[ERROR 4] - ' + 'Column name of the third column in MS1 peak table is not "rt"',
                    5: '[WARNING 5] - ' + 'The maximum value of "rt" is less than 60, the unit of rt may not be "second"',
                    6: '[WARNING 6] - ' + 'At least one MS1 peak in the table is found to have more than 50% of zero values among its peak abundances',
                    7: '[ERROR 7] - ' + 'Column name of the first column in Sample Information is not "sample.name"',
                    8: '[ERROR 8] - ' + 'Column name of the second column in Sample Information is not "group"',
                    9: '[ERROR 9] - ' + 'Sample names in MS1 peak table and Sample Information are different'};


@Component({
  selector: 'app-file-upload',
  templateUrl: './file-upload.component.html',
  styleUrls: ['./file-upload.component.css']
})
export class FileUploadComponent implements OnInit {
  private files: UploadFile[] = [];
  public uploadInput: EventEmitter<UploadInput>;
  public options: UploaderOptions;
  private humanizeBytes: Function = humanizeBytes;
  private fileId2Name: Object = {};
  public maxSize: number = 500;
  public disablePrNaInput: boolean = false;
  public projectName: string = '';
  public projectNameInfo: string = '';
  private currentUser: object = {};
  // private headers: Headers = new Headers();
  public totalSize: number = 0;
  @Output() ms1UploadDone = new EventEmitter<string>();
  private fileId2OutInfo: Object = {};
  private fileContentInfo: Object = {'peak_size': {'pos': 0, 'neg': 0},
                                     'sample_size': {'pos': 0, 'neg': 0}};
  
  public polarity: string = 'pos';
  public fileId2Pol: Object = {};
  public pol2FileId: Object = {'pos': [], 'neg': []};
  private fileId2Type: Object = {};
  public fileType2Id: Object = {'ms1_data': [], 'sample_info': [], 'ms2_data': []};
  public fileTypeMapping: Object = {'ms1_data': 'MS1 Peak Table', 'ms2_data': 'MS/MS Data',
                                    'sample_info': 'Sample Info'};
  private fileStatus2Id: Object = {'valid': [], 'warning': [], 'error': [], 'checking': []};
  private uploadedFileId: string[] = [];
  // private ms1FileId: string = '';
  private sampleInfoFileId: string = '';
  private sampleInfoUploadTimer: number = 0;
  public sampleInfoStatus: number[] = [];
  public ms1StoredUrl: string = '';
  // private peakProfileUrl: any;
  public hintUploadInfo = '';

  public ifHasCorrectUploadOrder = true;  

  private errorCode = ERROR_CODE;
  private errorCodeKey = Object.keys(this.errorCode).filter(ele => ele !== '0');


  @ViewChild('sampleInfoInput') sampleInfoInput: ElementRef;
  @ViewChild('ms1PosInput') ms1PosInput: ElementRef;
  @ViewChild('ms2PosInput') ms2PosInput: ElementRef;
  @ViewChild('ms1NegInput') ms1NegInput: ElementRef;
  @ViewChild('ms2NegInput') ms2NegInput: ElementRef;
  @ViewChild('fileInputForm') fileInputForm: ElementRef;
  
  constructor(
    private userService: UserService,
    private fileService: FileService,
    private santizer: DomSanitizer
  ) {
    this.options = { concurrency: 1 };
    this.uploadInput = new EventEmitter<UploadInput>();
   }

  ngOnInit() {
    this.currentUser = this.userService.getCurrentUser();
    // this.polarityHelp.nativeElement.classList.add('popover');
    $("#sampleInfo-help").popover({
      toggle: "popover",
      placement: "right",
      trigger: "hover",
      html: true,
      content: "<p>Upload one 'sample information' file for both positive and negative ionization data sets</p>"
    });
    $("#polarity-help").popover({
      toggle: "popover",
      placement: "right",
      trigger: "hover",
      html: true,
      content: "<p>Select the <b>polarity mode</b> according to your experimental condition. " +
                "If you want to process both positive and negative data sets together, " + 
                "upload positive data first, then switch to upload negative mode data. " +
                "All data files should be uploaded to the same project name.</p>"
    });
    $("#projectName-help").popover({
      toggle: "popover",
      placement: "right",
      trigger: "hover",
      html: true,
      content: "<p>Enter a name for the project. If you want to process both positive and negative data sets together, " +
               "all data files should be uploaded to the same project name.</p>"
    });
  }

  // for ngx-uploader
  onUploadOutput(output: UploadOutput): void {
    // console.log(output);
    this.files = this.fileService.onUploadOutput(output, this.files);
    if (output.type === 'done') {
      let currentFileId: string = output.file.id;
      this.uploadedFileId.push(currentFileId);
      let res = output.file.response;
      let error_code: number[] = JSON.parse(res['validity']);  // an array
      this.fileId2OutInfo[currentFileId] = {'file_type': res['file_type'], 'status': 'checking', 'error_code': error_code};
      if (res['file_type'] === 'ms1_data' || res['file_type'] === 'sample_info') {
        if (res['file_type'] === 'ms1_data') {
          let file_info = JSON.parse(res['file_info']);
          let _pol: string = this.fileId2Pol[currentFileId];
          this.fileContentInfo['sample_name_in_ms1'] = file_info['sample_name'];
          this.fileContentInfo['peak_size'][_pol] = file_info['peak_size'];
          this.fileContentInfo['sample_size'][_pol] = file_info['sample_size'];
          this.ms1StoredUrl = res['stored_url'];
          error_code = this.fileService.compareSampleName(error_code, currentFileId, this.fileContentInfo);
          this.ms1UploadDone.emit(_pol);
        }
        else if (res['file_type'] === 'sample_info') {
          // console.log(this.fileId2OutInfo);
          this.sampleInfoUploadTimer += 1;
          let file_info = JSON.parse(res['file_info']);
          let groups = Object.keys(file_info);
          this.fileContentInfo['groups'] = groups;
          let ms1Ids = this.fileType2Id['ms1_data'];
          if (this.sampleInfoUploadTimer >= 2) {
            this.updateExistMs1();
          }
        }
      }
      this.fileId2OutInfo[currentFileId]['error_code'] = error_code;
      this.fileId2OutInfo = this.errorCode2Status(currentFileId, this.fileId2OutInfo);
      this.fileStatus2Id = this.updateStatus2Id(this.fileId2OutInfo);
    }
  }

  startUpload(): void {
    let idOfPol = this.pol2FileId[this.polarity];
    let fileName2IdOfPol = this.fileService.getName2IdOfPol(this.fileId2Name, idOfPol);
    let fileId2TypeOfPol = this.fileService.getId2TypeOfPol(this.fileId2Type, idOfPol);
    let selectedTypeOfPol:string[] = Object.keys(fileId2TypeOfPol).map(ele => fileId2TypeOfPol[ele]);
    let hasMs1 = selectedTypeOfPol.includes('ms1_data');
    let hasMs2 = selectedTypeOfPol.includes('ms2_data');
    let pol = this.polarity;
    let polNameMapping = {'pos': 'positive', 'neg': 'negative'};
    if (hasMs1 && hasMs2) {
      this.hintUploadInfo = '';
      this.disablePrNaInput = true;
      this.fileService.startUpload(this.currentUser, fileId2TypeOfPol, this.polarity,
                                   fileName2IdOfPol, this.projectName, this.uploadInput);
    }
    else if (!hasMs1 && !hasMs2) {
      this.hintUploadInfo = 'No MS1 and MS/MS data for ' + polNameMapping[pol] + ' polarity mode.';
    }
    else if (!hasMs1) {
      this.hintUploadInfo = 'No MS1 data for ' + polNameMapping[pol] + ' polarity mode.';
    }
    else {
      this.hintUploadInfo = 'No MS/MS data for ' + polNameMapping[pol] + ' polarity mode.';
    }

  }

  cancelUpload(id: string): void {
    this.fileService.cancelUpload(id, this.uploadInput);
  }

  removeFile(id: string): void {
    let file:UploadFile;
    let fileType = this.fileId2Type[id];
    let pol = this.fileId2Pol[id];
    this.files.forEach(function(ele) {
      if (ele.id === id) {
        file = ele;
      }
    });
    if (file.responseStatus === 201) {
      let fileId = file.response.id;
      let token = this.currentUser['token'];
      this.fileService.deleteFile(token, fileId);
    }
    this.fileService.removeFile(id, this.uploadInput);
    this.fileType2Id[fileType] = this.fileType2Id[fileType].filter(ele => ele !== id);
    delete this.fileId2Type[id];
    if (this.uploadedFileId.includes(id)) {
      this.uploadedFileId = this.uploadedFileId.filter(ele => ele!==id);
    }

    if (this.fileId2OutInfo.hasOwnProperty(id)) {
      let status = this.fileId2OutInfo[id].status;
      if (this.fileStatus2Id[status].includes(id)) {
        this.fileStatus2Id[status] = this.fileStatus2Id[status].filter(ele => ele !== id);
      }
      delete this.fileId2OutInfo[id];
    }

    if (this.fileId2Pol.hasOwnProperty(id)) {

      if (this.pol2FileId[pol].includes(id)) {
        this.pol2FileId[pol] = this.pol2FileId[pol].filter(ele => ele !== id);
      }
      delete this.fileId2Pol[id];
    }

    if (fileType === 'sample_info') {
      this.sampleInfoFileId = '';
      this.sampleInfoStatus = [];
      this.sampleInfoInput.nativeElement.value = '';
    }
    else if (fileType === 'ms1_data') {
      if (pol === 'pos') {
        this.ms1PosInput.nativeElement.value = '';
      }
      else {
        this.ms1NegInput.nativeElement.value = '';
      }
    }
    else {
      let fileIdOfPol = this.pol2FileId[pol];
      let fileIdOfPolMs2 = fileIdOfPol.filter(ele => this.fileType2Id['ms2_data'].includes(ele));
      if (fileIdOfPolMs2.length === 0) {
        if (pol == 'pos') {
          this.ms2PosInput.nativeElement.value = '';
        }
        else {
          this.ms2NegInput.nativeElement.value = '';
        }
      }
    }

    setTimeout(() => {
      this.updateSelectedFileSize();
    }, 50);
  }

  removeAllFiles(): void {
    this.projectNameInfo = '';
    this.sampleInfoFileId = '';
    this.fileId2Type = {};
    this.fileType2Id = {'ms1_data': [], 'sample_info': [], 'ms2_data': []};
    this.sampleInfoUploadTimer = 0;
    this.sampleInfoStatus = [];
    this.fileId2Name = {};
    this.fileService.removeAllFiles(this.uploadInput);
    this.files = [];
    this.fileContentInfo = {'peak_size': {'pos': 0, 'neg': 0},
                            'sample_size': {'pos': 0, 'neg': 0}};
    this.fileId2OutInfo = {};
    this.uploadedFileId = [];
    this.fileStatus2Id = {'valid': [], 'warning': [], 'error': [], 'checking': []};
    this.fileId2Pol = {};
    this.pol2FileId = {'pos': [], 'neg': []};
    this.updateSelectedFileSize();
    this.fileInputForm.nativeElement.reset();
    this.polarity = '0';
    this.disablePrNaInput = false;
    setTimeout(() => {
      this.selectPolarity('pos');
    }, 5);
  }

  onChange(event: EventTarget, type: string, pol: string): void {
    setTimeout(() => {
      let eventObj: MSInputMethodContext = <MSInputMethodContext> event;
      let target: HTMLInputElement = <HTMLInputElement> eventObj.target;
      let fileNum: number = target.files.length;
      let currentFiles: UploadFile[] = this.files.slice(this.files.length - fileNum);
      let id2type: Object = {};
      let type2id: Object = {};
      let pol2id: Object = {};
      let id2pol: Object = {};
      let id2name: Object = {};
      let _sampleNames = [];
      let _sampleInfoStatus = [];
      currentFiles.forEach(function(item) {
        id2type[item.id] = type;
        id2pol[item.id] = pol;
        id2name[item.id] = item.name;
        if (type2id[type] === undefined) {
          type2id[type] = [];
        }
        if (pol2id[pol] === undefined) {
          pol2id[pol] = [];
        }
        type2id[type].push(item.id);
        pol2id[pol].push(item.id);

      });

      this.fileId2Name = {...this.fileId2Name, ...id2name};
      this.fileId2Type = {...this.fileId2Type, ...id2type};
      this.fileType2Id[type] = this.fileType2Id[type].concat(type2id[type]);
      this.fileId2Pol = {...this.fileId2Pol, ...id2pol};
      this.pol2FileId[pol] = this.pol2FileId[pol].concat(pol2id[pol])

      if (type === 'sample_info') {
        let currentFile = target.files[0];
        let sampleInfo = this.fileService.getSampleName(currentFile);
        sampleInfo.then(res => {
          if (res[0] !== 'sample.name') {
            this.sampleInfoStatus.push(7);
          }
          else {
            this.fileContentInfo['sample_name_in_sample_info'] = res.slice(1);
          }
        });
        this.sampleInfoFileId = currentFiles[0].id;
      }

      this.updateSelectedFileSize();
    }, 50)
  }

  updateSelectedFileSize() {
    let _totalSize: number = 0;
    this.files.forEach(function(ele) {
      _totalSize += ele.size;
    });
    this.totalSize = _totalSize;
  }

  errorCode2Status(fileId: string, fileId2OutInfo: Object) {
    let error_code: number[] = fileId2OutInfo[fileId]['error_code'].slice(0);  // an array
    let file_type: string = fileId2OutInfo[fileId]['file_type'];

    if (error_code.includes(0)) {
      if (error_code.length > 1) {
        fileId2OutInfo[fileId]['status'] = 'warning';
        fileId2OutInfo[fileId]['error_code'] = error_code.filter(ele => ele!==0);
      }
      else {
        // valid file
        fileId2OutInfo[fileId]['status'] = 'valid';
      }
    }
    else {
      fileId2OutInfo[fileId]['status'] = 'error';
    }
    return fileId2OutInfo;
  }

  updateStatus2Id(fileId2OutInfo: Object) {
    let _fileStatus2Id: Object = {'valid': [], 'warning': [], 'error': [], 'checking': []};
    Object.keys(fileId2OutInfo).forEach(id => {
      let _status = fileId2OutInfo[id]['status'];
      _fileStatus2Id[_status].push(id);
    });
    return _fileStatus2Id;
  }

  updateExistMs1() {
    let ms1Ids = this.fileType2Id['ms1_data'];
    if (ms1Ids.length !== 0) {
      let _com = this.fileService.compareSampleName;
      let _errorCode2Status = this.errorCode2Status;
      let _update = this.updateStatus2Id;
      let _uploadedId = this.uploadedFileId;
      ms1Ids.forEach((id, index) => {
        if (_uploadedId.includes(id)) {
          let error_code: number[] = this.fileId2OutInfo[id]['error_code'];
          this.fileId2OutInfo[id]['status'] = 'checking';
          error_code = _com(error_code, id, this.fileContentInfo);
          this.fileId2OutInfo[id]['error_code'] = error_code;
          this.fileId2OutInfo = _errorCode2Status(id, this.fileId2OutInfo);
          this.fileStatus2Id = _update(this.fileId2OutInfo);
        }
      })
    }
  }

  selectPolarity(polarityType: string) {
    this.polarity = polarityType;
    this.ifHasCorrectUploadOrder = true;
    let samInId = this.sampleInfoFileId;
    if (samInId) {
      if (polarityType === 'neg') {
        this.pol2FileId['pos'] = this.pol2FileId['pos'].filter(ele => ele!==samInId);
        if (this.pol2FileId['pos'].length !== 0 && this.files.length !== this.uploadedFileId.length) {
          this.ifHasCorrectUploadOrder = false;
        }
      }
      else {
        this.pol2FileId['neg'] = this.pol2FileId['neg'].filter(ele => ele!==samInId);
        if (this.pol2FileId['neg'].length !== 0 && this.files.length !== this.uploadedFileId.length) {
          this.ifHasCorrectUploadOrder = false;
        }
      }
      this.pol2FileId[polarityType].push(samInId);
      this.fileId2Pol[samInId] = polarityType;
    }
  }
}
