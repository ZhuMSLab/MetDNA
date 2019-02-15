// created by Xin Xiong<onlybelter@outlook.com>, https://github.com/OnlyBelter

import { Injectable, EventEmitter } from '@angular/core';
import { Headers, Http, RequestOptions, ResponseContentType } from "@angular/http";
import { Observable } from 'rxjs';
import { UploadOutput, UploadInput, UploadFile, humanizeBytes, UploaderOptions} from 'ngx-uploader';
import { ProjectService } from './project.service';
import { fileApi } from '../_data_model/index';

// const URL = 'http://192.168.201.211:8024/files/';
// const URL = djangoApi
// const mediaApi = 'http://localhost:8000/media/';

@Injectable()
export class FileService {
  private fileUrl: string = fileApi;

  constructor(
    private projectService: ProjectService,
    private http: Http
  ) { }

  onUploadOutput(output: UploadOutput, files: UploadFile[]): UploadFile[] {
    if (output.type === 'allAddedToQueue') { // when all files added in queue
      // uncomment this if you want to auto upload files when added
      // const event: UploadInput = {
      //   type: 'uploadAll',
      //   url: '/upload',
      //   method: 'POST',
      //   data: { foo: 'bar' },
      //   concurrency: 0
      // };
      // this.uploadInput.emit(event);
    } else if (output.type === 'addedToQueue'  && typeof output.file !== 'undefined') { // add file to array when added
      files.push(output.file);
    } else if (output.type === 'uploading' && typeof output.file !== 'undefined') {
      const index = files.findIndex(file => typeof output.file !== 'undefined' && file.id === output.file.id);
      files[index] = output.file;
    } else if (output.type === 'removed') {
      files = files.filter((file: UploadFile) => file !== output.file);
    } else if (output.type === 'removedAll') {
      files = [];
    }
    return files;
  }

  startUpload(user: Object, fileId2Type: Object, currentPol: string, fileName2Id: Object, 
              projectName: string, uploadInput: EventEmitter<UploadInput>): void {
    let token = user['token'];
    let checkProjectNames = this.projectService.getProjectByName(token, projectName);
    // console.log(checkProjectNames);
    checkProjectNames.then(res => {
      if (res[0]) {
        this.uploadEmitter(token, this.fileUrl, fileId2Type, currentPol, fileName2Id, projectName, uploadInput);
      }
      else {
        // project name not exist
        // console.log('I am null');
        if (token && projectName) {
          let projectInfo = this.createProject(token, projectName);
          projectInfo.then(res => {
            if (res.status === 201) {
              // create new project successfully
              this.projectService.saveNewProjectName(projectName);
              this.uploadEmitter(token, this.fileUrl, fileId2Type, currentPol, fileName2Id, projectName, uploadInput);
              // this.projectName = '';
              // this.projectNameInfo = 'This project has created successfully,' +
              //                          ' press "Reset" to start a new project!';
            }
          })
        }
      }
    })
  }

  uploadEmitter(token: string, url: string, fileId2Type: Object, currentPol: string,
    fileName2Id: Object, projectName: string, uploadInput: EventEmitter<UploadInput>) {
    let id2Type = JSON.stringify(fileId2Type);
    let name2Id = JSON.stringify(fileName2Id);
    // let projectName = projectName;
    const event: UploadInput = {
    type: 'uploadAll',
    // url: 'http://ngx-uploader.com/upload',
    url: url,
    method: 'POST',
    headers: { 'Authorization': 'JWT ' + token },
    data: { fileId2Type: id2Type, currentPol: currentPol,
            fileName2Id: name2Id, projectName: projectName },
    // withCredentials: true
    };
    uploadInput.emit(event);
  }

  cancelUpload(id: string, uploadInput: EventEmitter<UploadInput>): void {
    uploadInput.emit({ type: 'cancel', id: id });
  }

  removeFile(id: string, uploadInput: EventEmitter<UploadInput>): void {
    uploadInput.emit({ type: 'remove', id: id });
  }

  removeAllFiles(uploadInput: EventEmitter<UploadInput>): void {
    
    uploadInput.emit({ type: 'removeAll' });
  }

  deleteFile(token: string, fileId: number) {
    const url = this.fileUrl + fileId + '/';
    let headers = new Headers();
    headers.append('Authorization', 'JWT ' + token);
    // headers.append('Content-Type', 'application/json');
    let options = new RequestOptions({ headers: headers });
    return this.http.delete(url, options)
                    .toPromise()
                    .then(() => null)
                    .catch(this.handleError);
  }

  createProject(token: string, projectName: string) {
    return this.projectService.createProject(token, projectName);
  }

  getPeakProfile(storedUrl: string): Promise<File> {
    let file: File;
    let headers = new Headers({'Content-Type': 'image/jpg'});
    return this.http.get(storedUrl, {headers: headers, responseType: ResponseContentType.Blob})
                    .toPromise()
                    .then(res => {
                        if (res.status===200) {
                        return res.blob() as File;
                      }
                        else { 
                          return file; 
                        }
                      })
                    .catch(this.handleError)
  }

  compareSampleName(errorCode: number[], ms1Id: string, fileContentInfo: Object) {
    let sample_name_in_ms1: string[] = fileContentInfo['sample_name_in_ms1'];
    let sample_name_in_sample_info: string[] = fileContentInfo['sample_name_in_sample_info'];
    let id = ms1Id;
    let _error_code: number[] = errorCode;
    let _error_code_temp = _error_code.slice(0);  // deep copy array
    if (sample_name_in_ms1 && sample_name_in_sample_info) {
      let str1 = JSON.stringify(sample_name_in_ms1.sort());
      let str2 = JSON.stringify(sample_name_in_sample_info.sort());
      if (str1 !== str2) {
        if (_error_code.includes(0)) {
          _error_code_temp = _error_code.filter(ele => ele !== 0).concat([9]);
        }
        else if (!_error_code.includes(9)) {
          _error_code_temp.push(9);
        }
      }
      else {
        if (_error_code.includes(9)) {
          _error_code_temp = _error_code.filter(ele => ele !== 9);
        }
        if ((_error_code_temp.length===1 && [5, 6].includes(_error_code_temp[0])) ||
          _error_code_temp.length===0) {
          _error_code_temp = _error_code_temp.concat([0]);
        }
      }
    }
    return _error_code_temp;
  }

  getSampleName(sampleInfoFile: File): Promise<string[]> {
    return new Promise((resolve, reject) => {
      let fr = new FileReader();
      fr.onloadend = () => {
        let rawData = fr.result.replace(/\r\n?/g, '\n').split('\n');
        let sampleName = [];
        rawData = rawData.map((x) => x.split(','));
        rawData.forEach((ele) => {
          if (ele.length === 2) {
            let name = ele[0].replace(/^\"|\"$/g, "");
            sampleName.push(name);
          }
        });
        resolve(sampleName);
      }
      fr.readAsText(sampleInfoFile);
    });
  }

  getName2IdOfPol(id2name: Object, idOfPol: string[]): Object {
    let result:Object[] = Object.keys(id2name)
                              .filter(key => idOfPol.includes(key))
                              .map(key => {
                                let obj = {};
                                obj[id2name[key]] = key;
                                return obj;
                              });
    return Object.assign({}, ...result);
  }

  getId2TypeOfPol(id2type: Object, idOfPol: string[]): Object {
    let result:Object[] = Object.keys(id2type)
                                 .filter(key => idOfPol.includes(key))
                                 .map(key => {
                                   let obj = {};
                                   obj[key] = id2type[key];
                                   return obj;
                                 });
    return Object.assign({}, ...result);
  }

  private handleError(error: any): Promise<any> {
    console.error('An error occurred', error);
    return Observable.of(error).toPromise();
  }

}
