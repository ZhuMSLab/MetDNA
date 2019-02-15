import { Injectable } from '@angular/core';
import { Headers, Http, RequestOptions } from "@angular/http";
import { paraApi } from '../_data_model/index';

import 'rxjs/add/operator/toPromise';


@Injectable()
export class ParameterService {
  private paraApi = paraApi;

  constructor(
    private http: Http,
  ) { }

  createOptions(token: string) {
    let headers = new Headers();
    headers.append('Authorization', 'JWT ' + token);
    let options = new RequestOptions({ headers: headers });
    return options;
  }

  createProjectQueue(token: string, projectName: string, paras: Object):Promise<number> {
    const url = this.paraApi;
    let options = this.createOptions(token);
    options.headers.append('Content-Type', 'application/json');
    return this.http.
           post(url, JSON.stringify({project_name: projectName, paras: paras}), options)
           .toPromise()
           .then(res => {
            //  console.log(res);
             return res.status;
           })
           .catch(this.handleError);
  }

  private handleError(error: any): Promise<number> {
    console.error('An error occurred', error); // for demo purposes only
    return Promise.reject(error.status);
  }




}
