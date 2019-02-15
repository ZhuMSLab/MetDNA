import { Component, OnInit } from '@angular/core';
import { UserService } from './_service/index';
import { frontPort } from './_data_model/index';

@Component({
  selector: 'app-root',
  templateUrl: './app.component.html',
  styleUrls: ['./app.component.css']
})
export class AppComponent implements OnInit{
  title = 'MetDNA';
  subtitle = 'Metabolite annotation and Dysregulated Network Analysis';
  private LOGO = '/assets/logo-mini.png';
  public frontPort = frontPort;

  constructor(private userService: UserService) { }

  ngOnInit() {
    let currentUser = this.userService.getCurrentUser();
    let token: string = '';
    if (currentUser) {
      token = currentUser['token'];
      let checkTokenInfo = this.userService.checkUserToken(token);
      // console.log(checkTokenInfo);
      checkTokenInfo.then(res => {
        if (res===400) {
          this.userService.userLogout();
        }
      })
    }

  }

  
  // @HostListener('window:beforeunload')
  // private onUnload(): void {
  //   // this.userService.userLogout();
  //   alert('are you sure?');
  // }
  
}
