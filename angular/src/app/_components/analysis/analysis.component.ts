// created by Xin Xiong<onlybelter@outlook.com>, https://github.com/OnlyBelter -->

import { Component, OnInit, EventEmitter, ViewChild  } from '@angular/core';
import {FormGroup, FormBuilder, Validators} from '@angular/forms';
import {DomSanitizer} from '@angular/platform-browser';
import { FileUploadComponent } from '../../_components/file-upload/file-upload.component';
import { UserService, FileService } from '../../_service/index';
import { frontPort } from '../../_data_model/index';

// get MS1 peak intensity profile figure in this component
const peakIntFileName = 'peak_intensity_profile.png';
@Component({
  selector: 'app-analysis',
  templateUrl: './analysis.component.html',
  styleUrls: ['./analysis.component.css']
})

export class AnalysisComponent implements OnInit {

  private email: string = '';
  private hideEmailForm: boolean = false;
  private disableWindow: Object = {1: false, 2: true, 3: true, 4: true};
  private currentUser: object = {};
  private window: number = 1;
  private finalSubmitted: boolean = false;
  private frontPort: string = frontPort;
  // private showSubmittedInfo: boolean = false;
  public peakProfileUrlPos: any;
  public peakProfileUrlNeg: any;
  public polarityMode: string = 'pos';  // pos/neg/both
  public ms2FileCount: Object = {'pos': 0, 'neg': 0};
  // public maxSize: number = 500;  // the max size of total upload files is 500Mb
  // this way can call other component's function directly in this component
  @ViewChild(FileUploadComponent)
  private fileComponent: FileUploadComponent;

  constructor(
    private fb: FormBuilder,
    private userService: UserService,
    private fileService: FileService,
    private santizer: DomSanitizer
  ) { }

  ngOnInit() {

    let user = this.userService.getCurrentUser();
    if (user) {
      this.currentUser = user;
      this.email = this.currentUser['email'];
      let token = this.currentUser['token'];
      let checkUserInfo = this.userService.checkUserToken(token);
      checkUserInfo.then(res => {
        if (res===400) {
          this.userService.userLogout();
        }
        else if (res===200) {
          this.hideEmailForm = true;
        }
      })
    }
  }

  // getPeakProfileUrl(pol: string) { //
  //   // let User = this.userService.getUserByEmail(this.email);

  // }

  getPolarityMode(pol2id: Object) {
    let polMode: string = '';
    let posFileCount = pol2id['pos'].length;
    let negFileCount = pol2id['neg'].length;
    if (posFileCount >= 2 && negFileCount >= 2) {
      // at least two files in this polarity mode
      polMode = 'both';
    }
    else if (posFileCount >= 2) {
      polMode = 'pos';
    }
    else {
      polMode = 'neg';
    }
    return polMode;
  }

  changeWindow(id: number) {
    if (this.disableWindow[id] === false) {
      this.window = id;
    }
  }

  getMs2FileCount(ms2Ids: string[], pol2id: Object) {
    let _ms2FileCount: Object = {'pos': 0, 'neg': 0};
    Object.keys(pol2id).forEach((ele) => {
      let _IdOfPol = pol2id[ele];
      let _ms2IdOfPol = _IdOfPol.filter(id => ms2Ids.includes(id));
      _ms2FileCount[ele] = _ms2IdOfPol.length;
    })
    return _ms2FileCount;
  }

  openNextWindow(id: number) {
    if (id === 3) {
      let _pol2id = this.fileComponent.pol2FileId;
      let _ms2Ids = this.fileComponent.fileType2Id['ms2_data'];
      this.polarityMode = this.getPolarityMode(_pol2id);
      // console.log(this.polarityMode);
      this.ms2FileCount = this.getMs2FileCount(_ms2Ids, _pol2id);
      // console.log(this.ms2FileCount);
    }

    this.window = id;
    this.disableWindow[id] = false;
  }

  // an EventEmitter defined in analysis-para-setting component
  finalSubmit(submit: boolean) {
    if (submit) {
      this.window = 2;
      // this.disableWindow[2] = true;
      this.disableWindow[3] = true;
      this.disableWindow[4] = true;
      this.fileComponent.removeAllFiles();
      // console.log('after final submit...');
    }
  }
  
  // an EventEmitter defined in file-upload component
  // EventEmitter can pass data from child component to parent component
  ms1UploadDone(pol: string) {
    let ms1Url: string = this.fileComponent.ms1StoredUrl;
    pol = pol.toUpperCase();
    // console.log(pol);
    if (ms1Url) {
      let _tem: string[] = ms1Url.split('/');
      _tem = _tem.slice(0, _tem.length-2);
      // console.log(_tem);
      let _newPeakProfileName = peakIntFileName.replace('.png', `_${pol}.png`);
      _tem.push(_newPeakProfileName);
      let baseUrl = _tem.join('/');
      // change port from back-end port to frontend port to avoid Access-Control-Allow-Origin
      baseUrl = baseUrl.replace(/\d{4}/, this.frontPort);  
      // console.log(baseUrl);
      let _peakProfile = this.fileService.getPeakProfile(baseUrl);
      // console.log(_peakProfile);
      _peakProfile.then(res => {
        // console.log(res);
        let urlCreator = window.URL;
        if (pol === 'POS') {
          this.peakProfileUrlPos = this.santizer.bypassSecurityTrustUrl(urlCreator.createObjectURL(res));
        }
        else {
          this.peakProfileUrlNeg = this.santizer.bypassSecurityTrustUrl(urlCreator.createObjectURL(res));
        }
        // console.log(this.peakProfileUrl);
      })
    }
  }

}


