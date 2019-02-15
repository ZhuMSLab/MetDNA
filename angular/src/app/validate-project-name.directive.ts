import { Directive, Input } from '@angular/core';
import { NG_ASYNC_VALIDATORS, Validator, AbstractControl } from '@angular/forms';
import { UserService, ProjectService } from './_service/index';


@Directive({
  selector: '[validateProjectName]',
  providers: [{provide: NG_ASYNC_VALIDATORS, useExisting: ValidateProjectNameDirective, multi: true}]
})
export class ValidateProjectNameDirective implements Validator{
  // @Input() validateProjectName: string;
  private projectNames: Promise<string[]>;

  constructor(
    private userService: UserService,
    private projectService: ProjectService,
  ) { }

  // https://stackoverflow.com/a/36238205/2803344
  validate(control: AbstractControl) {
    return new Promise(resovle => {
      // console.log(control);
      // console.log(this.projectNames);

      if (this.projectNames === undefined) {
        let user = this.userService.getCurrentUser();
        let token = user['token'];
        // get all current projects' name of this user
        this.projectNames = this.projectService.getUserProjectNames(token);
      }
      let isValid = true;
      this.projectNames.then(res => {
        // console.log(res);
        let newProjectName = this.projectService.getNewProjectName();
        let totalProjectName: string[] = res.slice();
        if (newProjectName) {
          totalProjectName = totalProjectName.concat(newProjectName);
        }
        // console.log(totalProjectName);
        if (totalProjectName.includes(control.value)) {
          isValid = false;
        }
        // console.log(isValid);
        if (isValid) {
          resovle(null);
        }
        else {
          resovle({projectName: {valid: false}});
        }
      });

    })
  }

}
