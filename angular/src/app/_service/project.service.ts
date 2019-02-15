// created by Xin Xiong<onlybelter@outlook.com>, https://github.com/OnlyBelter
import { Injectable } from '@angular/core';
import { Headers, Http, RequestOptions, Response } from "@angular/http";
import { Project } from '../_data_model/index';
import { projectApi } from '../_data_model/index';

import 'rxjs/add/operator/toPromise';

@Injectable()
export class ProjectService {
  private projectApi = projectApi;

  constructor(
    private http: Http,
  ) { }

  createOptions(token: string) {
    let headers = new Headers();
    headers.append('Authorization', 'JWT ' + token);
    let options = new RequestOptions({ headers: headers });
    return options;
  }

  getProjectByName(token: string, projectName: string): Promise<Project> {
    const url = this.projectApi;
    let headers = new Headers();
    headers.append('Authorization', 'JWT ' + token);
    headers.append('Content-Type', 'application/json');
    let options = new RequestOptions({ headers: headers, params: {'project_name': projectName} });
    return this.http.get(url, options)
               .toPromise()
               .then(res => {
                 return res.json().results as Project;
              })
               .catch(this.handleError);
  }

  getProject(token: string): Promise<Project[]> {
    const url = this.projectApi;
    // console.log('here is in user.service.ts', url);
    let options = this.createOptions(token);
    return this.http.get(url, options)
               .toPromise()
               .then(res => {
                 return res.json().results as Project[];
              })
               .catch(this.handleError);
  }

  getUserProjectNames(token: string): Promise<string[]> {
    // get this user's total projects' name
    localStorage.removeItem('newProjectName');
    let projectInfo = this.getProject(token);
    // console.log(projectInfo);
    let projectNames = [];
    return projectInfo.then(res => {
      // console.log(res);
      res.forEach(function(ele) {
        projectNames.push(ele['project_name']);
      });
      // console.log(projectNames);
      return projectNames;
    });
  }

  createProject(token: string, projectName: string):Promise<Response> {
    const url = this.projectApi;
    let options = this.createOptions(token);
    options.headers.append('Content-Type', 'application/json');
    let projectHash = this.generateProjectId();
    return this.http.
           post(url, JSON.stringify({project_hash: projectHash, 
                                project_name: projectName}), options)
           .toPromise()
           .then(res => res)
           .catch(this.handleError);
  }

  saveNewProjectName(name: string) {
    let newProjectName: string[] = [name];
    let prevailNewProjectName: string[] = JSON.parse(localStorage.getItem('newProjectName'));
    if (prevailNewProjectName) {
      newProjectName = [name].concat(prevailNewProjectName);
    }
    localStorage.setItem('newProjectName', JSON.stringify(newProjectName));
  }

  getNewProjectName(): string[] {
    return JSON.parse(localStorage.getItem('newProjectName'));
  }


  private handleError(error: any): Promise<any> {
    console.error('An error occurred', error); // for demo purposes only
    return Promise.reject(error.message || error);
  }

    // https://stackoverflow.com/a/2117523/2803344
    generateProjectId = function() {
      return 'xxxx8xyx'.replace(/[xy]/g, function(c) {
        var r = Math.random() * 16 | 0, v = c == 'x' ? r : (r & 0x3 | 0x8);
        return v.toString(16);
      });
    }



}
